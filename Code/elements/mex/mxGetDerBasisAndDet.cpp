//----------------------------------------------------------------------------------------
// computeJacobianBasis_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [N, detJ] = computeJacobianBasis(J1, coords, weights)
//     J1      : 3 x Nnodes x nGauss
//     coords  : Nnodes x 3
//     weights : nGauss x 1  (or 1 x nGauss)
//     N       : 3 x Nnodes x nGauss
//     detJ    : nGauss x 1
//
// Build command:
//   mex -R2018a computeJacobianBasis_wrap.cpp
//
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgIdAndTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() in the gateway (MexFunction) and
//         with a C++ throw in the math kernel — see [NEW-1] below.
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         All input arrays copied into std::vector<double>; .data()
//         passed to the computational kernels.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         Strings use factory.createCharArray() instead.
//
// [NEW-1] inv3x3() called mexErrMsgIdAndTxt() directly — this is
//         illegal in the C++ MEX API (not declared) and also unsafe
//         because longjmp would skip C++ destructors including
//         std::vector. Fixed by:
//           - Removing the mexErrMsgIdAndTxt call from inv3x3()
//           - Having inv3x3() throw std::runtime_error instead
//           - Catching it in operator()() and re-throwing via throwError()
//         This keeps the math kernels free of any MATLAB API dependency
//         and ensures all destructors (RAII) run correctly.
//
// [NEW-2] mxGetPr() / mxCreateNumericArray() replaced with
//         TypedArray<double> copies (input) and
//         factory.createArray<double>({...}) (output).
//         3D arrays use the initializer_list {d0, d1, d2} overload.
//
// [NEW-3] mxGetDimensions() / mxGetNumberOfDimensions() replaced with
//         Array::getDimensions() which returns a std::vector<size_t>.
//
// [NEW-4] Output written into std::vector<double> buffers during the
//         Gauss loop, then copied into TypedArray at the end.
//         This avoids the proxy-element issue for the 3D output array N.
//
// [NEW-5] std::size_t / mwSize used for all array sizes and loop bounds.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>   // std::runtime_error — used by inv3x3 [NEW-1]

using namespace matlab::data;
using matlab::mex::ArgumentList;

//========================================================================================
// Pure math kernels — NO MATLAB API dependency
// [NEW-1] Errors are reported via C++ exceptions, not mexErrMsgIdAndTxt
//========================================================================================

// Determinant of column-major 3×3 matrix
static double det3x3(const double* M)
{
    return M[0] * (M[4]*M[8] - M[5]*M[7])
         - M[1] * (M[3]*M[8] - M[5]*M[6])
         + M[2] * (M[3]*M[7] - M[4]*M[6]);
}

// Inverse of column-major 3×3 matrix
// [NEW-1] Throws std::runtime_error instead of calling mexErrMsgIdAndTxt —
//         keeps the kernel MATLAB-free and lets RAII destructors unwind.
static void inv3x3(const double* M, double* invM)
{
    const double det = det3x3(M);
    if (std::fabs(det) < 1e-15)
        throw std::runtime_error("Jacobian matrix is singular.");

    const double invDet = 1.0 / det;

    invM[0] =  (M[4]*M[8] - M[5]*M[7]) * invDet;
    invM[1] = -(M[1]*M[8] - M[2]*M[7]) * invDet;
    invM[2] =  (M[1]*M[5] - M[2]*M[4]) * invDet;
    invM[3] = -(M[3]*M[8] - M[5]*M[6]) * invDet;
    invM[4] =  (M[0]*M[8] - M[2]*M[6]) * invDet;
    invM[5] = -(M[0]*M[5] - M[2]*M[3]) * invDet;
    invM[6] =  (M[3]*M[7] - M[4]*M[6]) * invDet;
    invM[7] = -(M[0]*M[7] - M[1]*M[6]) * invDet;
    invM[8] =  (M[0]*M[4] - M[1]*M[3]) * invDet;
}

// Multiply (3×N) × (N×3) → 3×3 (all column-major)
static void mul3xN_Nx3(const double* A, const double* B, double* out, int N)
{
    for (int col = 0; col < 3; ++col)
        for (int row = 0; row < 3; ++row) {
            double sum = 0.0;
            for (int k = 0; k < N; ++k)
                sum += A[row + 3*k] * B[k + N*col];
            out[row + 3*col] = sum;
        }
}

// Multiply (3×3) × (3×N) → 3×N (all column-major)
static void mul3x3_3xN(const double* A, const double* B, double* out, int N)
{
    for (int col = 0; col < N; ++col)
        for (int row = 0; row < 3; ++row) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k)
                sum += A[row + 3*k] * B[k + 3*col];
            out[row + 3*col] = sum;
        }
}

//========================================================================================
// MexFunction
//========================================================================================
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-3] Read input arrays into std::vector via TypedArray copy
        // -----------------------------------------------------------------------
        const TypedArray<double> J1_arr      = inputs[0];
        const TypedArray<double> coords_arr  = inputs[1];
        const TypedArray<double> weights_arr = inputs[2];

        // [NEW-3] getDimensions() replaces mxGetDimensions / mxGetNumberOfDimensions
        const auto dimsJ1 = J1_arr.getDimensions();
        const int  Nnodes = static_cast<int>(dimsJ1[1]);
        const int  nGauss = (dimsJ1.size() >= 3) ? static_cast<int>(dimsJ1[2]) : 1;

        // [FIX-B] Copy into contiguous vectors — TypedArray proxies cannot be
        //         addressed directly for raw pointer arithmetic in the kernels
        std::vector<double> J1     (J1_arr.begin(),      J1_arr.end());
        std::vector<double> coords (coords_arr.begin(),  coords_arr.end());
        std::vector<double> weights(weights_arr.begin(), weights_arr.end());

        // -----------------------------------------------------------------------
        // [NEW-4] Allocate output buffers as std::vector — written by the loop,
        //         then copied into TypedArray after (avoids proxy write issues)
        // -----------------------------------------------------------------------
        const std::size_t szN    = static_cast<std::size_t>(3) *
                                   static_cast<std::size_t>(Nnodes) *
                                   static_cast<std::size_t>(nGauss);
        const std::size_t szDetJ = static_cast<std::size_t>(nGauss);

        std::vector<double> N_vec   (szN,    0.0);
        std::vector<double> detJ_vec(szDetJ, 0.0);

        // Temporary per-Gauss-point buffers
        double             J[9], invJ[9];
        std::vector<double> Ntemp(static_cast<std::size_t>(3 * Nnodes));

        // -----------------------------------------------------------------------
        // Main Gauss loop
        // [NEW-1] inv3x3 now throws std::runtime_error for singular Jacobian;
        //         caught here and re-thrown as a MATLAB exception so that all
        //         C++ destructors (std::vector etc.) unwind correctly.
        // -----------------------------------------------------------------------
        try {
            for (int gp = 0; gp < nGauss; ++gp) {
                const double* J1gp =
                    J1.data() + static_cast<std::size_t>(gp) *
                                static_cast<std::size_t>(3 * Nnodes);

                // J = J1(:,:,gp) * coords   (3×Nnodes × Nnodes×3 → 3×3)
                mul3xN_Nx3(J1gp, coords.data(), J, Nnodes);

                // detJ(gp) = det(J) * weight(gp)
                detJ_vec[static_cast<std::size_t>(gp)] =
                    det3x3(J) * weights[static_cast<std::size_t>(gp)];

                // invJ = J^{-1}  — may throw if singular
                inv3x3(J, invJ);

                // N(:,:,gp) = invJ * J1(:,:,gp)   (3×3 × 3×Nnodes → 3×Nnodes)
                mul3x3_3xN(invJ, J1gp, Ntemp.data(), Nnodes);

                // Store into flat output buffer
                const std::size_t offset =
                    static_cast<std::size_t>(gp) *
                    static_cast<std::size_t>(3 * Nnodes);
                for (int i = 0; i < 3 * Nnodes; ++i)
                    N_vec[offset + static_cast<std::size_t>(i)] = Ntemp[static_cast<std::size_t>(i)];
            }
        } catch (const std::runtime_error& e) {
            throwError("mexFunction:singularMatrix", e.what());
        }

        // -----------------------------------------------------------------------
        // [NEW-2] Build TypedArray outputs and copy results in
        // -----------------------------------------------------------------------

        // N : 3 × Nnodes × nGauss
        TypedArray<double> N_out =
            factory.createArray<double>(
                {3,
                 static_cast<std::size_t>(Nnodes),
                 static_cast<std::size_t>(nGauss)});
        std::copy(N_vec.begin(), N_vec.end(), N_out.begin());

        // detJ : nGauss × 1
        TypedArray<double> detJ_out =
            factory.createArray<double>(
                {static_cast<std::size_t>(nGauss), 1});
        std::copy(detJ_vec.begin(), detJ_vec.end(), detJ_out.begin());

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(N_out);
        outputs[1] = std::move(detJ_out);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 3)
            throwError("mexFunction:invalidNumInputs",
                       "3 inputs required: J1 (3 x Nnodes x nGauss), "
                       "coords (Nnodes x 3), weights (nGauss x 1).");

        if (outputs.size() > 2)
            throwError("mexFunction:invalidNumOutputs",
                       "2 outputs required: N (3 x Nnodes x nGauss), detJ (nGauss x 1).");

        // All three inputs must be real double arrays
        for (std::size_t i = 0; i < 3; ++i)
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("mexFunction:invalidInput",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a real double array.");

        // J1: first dim must be 3, at least 2D
        {
            const auto d = inputs[0].getDimensions();
            if (d.size() < 2 || d[0] != 3)
                throwError("mexFunction:invalidInput",
                           "J1 must be of size 3 x Nnodes x nGauss.");
        }

        // coords: must be 2D with 3 columns and Nnodes rows
        {
            const auto dJ1     = inputs[0].getDimensions();
            const auto dCoords = inputs[1].getDimensions();
            const std::size_t Nnodes = dJ1[1];
            if (dCoords.size() != 2 || dCoords[0] != Nnodes || dCoords[1] != 3)
                throwError("mexFunction:invalidInput",
                           "coords must be Nnodes x 3.");
        }

        // weights: must be a vector of length nGauss
        {
            const auto dJ1 = inputs[0].getDimensions();
            const std::size_t nGauss =
                (dJ1.size() >= 3) ? dJ1[2] : 1;
            const auto dW = inputs[2].getDimensions();
            const bool ok = (dW.size() >= 1) &&
                            ((dW[0] == nGauss && (dW.size() < 2 || dW[1] == 1)) ||
                             (dW.size() >= 2 && dW[1] == nGauss && dW[0] == 1));
            if (!ok)
                throwError("mexFunction:invalidInput",
                           "weights must be a vector of length nGauss.");
        }
    }

    // [FIX-E] createCharArray for strings; routes through MATLAB error()
    void throwError(const std::string& id, const std::string& msg)
    {
        getEngine()->feval(u"error", 0,
            std::vector<Array>{
                factory.createCharArray(id),
                factory.createCharArray(msg)
            });
    }
};

//----------------------------------------------------------------------------------------
