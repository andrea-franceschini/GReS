//----------------------------------------------------------------------------------------
// computeQuad9Derivatives_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   dN = computeQuad9Derivatives(coords)   % coords: n×2, dN: 2×9×n
//
// Build command:
//   mex -R2018a computeQuad9Derivatives_wrap.cpp
//
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
//         Generic IDs "Quad9Deriv:inputError" / "Quad9Deriv:outputError"
//         assigned since the original used the no-ID mexErrMsgTxt.
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         coords copied into std::vector<double>; kernel writes into a
//         std::vector<double> output buffer, then std::copy into TypedArray.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         All string arguments use factory.createCharArray() instead.
//
// [FIX-F] mwSize is declared in mex.h — NOT included by mex.hpp.
//         Every mwSize replaced with std::size_t (same underlying type
//         on all 64-bit platforms). #include <cstddef> provides it.
//
// [NEW-1] mxGetPr() / mxCreateNumericArray(3, {2,9,np}) replaced with
//         factory.createArray<double>({2, 9, np}) and std::copy.
//
// [NEW-2] mxGetM() / mxGetN() replaced with getDimensions()[0] / [1].
//
// [NEW-3] The pure computational function computeDerivatives() is kept
//         completely unchanged in logic — only the mwSize → std::size_t
//         substitution is applied. No MATLAB API dependency in the kernel.
//
// [NEW-4] 3D output indexing: dN[row + 2*col + 18*i]
//         The factory creates the array with dimensions {2, 9, np};
//         MATLAB stores 3D arrays in column-major order matching this
//         flat index exactly. Semantics identical to the original.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cstddef>   // std::size_t  [FIX-F]
#include <vector>
#include <string>
#include <algorithm> // std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// [NEW-3] Pure computational kernel — no MATLAB API dependency.
// [FIX-F] mwSize → std::size_t
//
// coords layout: column-major, np rows × 2 cols
//   xi  = coords[i]       (col 0)
//   eta = coords[i + np]  (col 1)
//
// dN layout: column-major 3D, dims {2, 9, np}
//   dN[row + 2*col + 18*i]
//     row = 0 → d/dxi,  row = 1 → d/deta
//     col = 0..8 (node index)
//     i   = 0..np-1 (point index)
//----------------------------------------------------------------------------------------
static void computeDerivatives(double* dN, const double* coords, std::size_t np)
{
    for (std::size_t i = 0; i < np; ++i) {
        const double xi  = coords[i];         // col 0
        const double eta = coords[i + np];    // col 1

        // 1D quadratic basis functions and their derivatives
        const double b1_xi  =  0.5 * xi  * (xi  - 1.0);
        const double b2_xi  =  1.0 - xi  * xi;
        const double b3_xi  =  0.5 * xi  * (xi  + 1.0);
        const double gb1_xi =  0.5 * (2.0 * xi  - 1.0);
        const double gb2_xi = -2.0 * xi;
        const double gb3_xi =  0.5 * (2.0 * xi  + 1.0);

        const double b1_eta  =  0.5 * eta * (eta - 1.0);
        const double b2_eta  =  1.0 - eta * eta;
        const double b3_eta  =  0.5 * eta * (eta + 1.0);
        const double gb1_eta =  0.5 * (2.0 * eta - 1.0);
        const double gb2_eta = -2.0 * eta;
        const double gb3_eta =  0.5 * (2.0 * eta + 1.0);

        const std::size_t base = 18 * i;   // stride per point

        // Node 1
        dN[base + 0 + 2 * 0] = gb1_xi * b1_eta;
        dN[base + 1 + 2 * 0] = b1_xi  * gb1_eta;
        // Node 2
        dN[base + 0 + 2 * 1] = gb3_xi * b1_eta;
        dN[base + 1 + 2 * 1] = b3_xi  * gb1_eta;
        // Node 3
        dN[base + 0 + 2 * 2] = gb3_xi * b3_eta;
        dN[base + 1 + 2 * 2] = b3_xi  * gb3_eta;
        // Node 4
        dN[base + 0 + 2 * 3] = gb1_xi * b3_eta;
        dN[base + 1 + 2 * 3] = b1_xi  * gb3_eta;
        // Node 5
        dN[base + 0 + 2 * 4] = gb2_xi * b1_eta;
        dN[base + 1 + 2 * 4] = b2_xi  * gb1_eta;
        // Node 6
        dN[base + 0 + 2 * 5] = gb3_xi * b2_eta;
        dN[base + 1 + 2 * 5] = b3_xi  * gb2_eta;
        // Node 7
        dN[base + 0 + 2 * 6] = gb2_xi * b3_eta;
        dN[base + 1 + 2 * 6] = b2_xi  * gb3_eta;
        // Node 8
        dN[base + 0 + 2 * 7] = gb1_xi * b2_eta;
        dN[base + 1 + 2 * 7] = b1_xi  * gb2_eta;
        // Node 9
        dN[base + 0 + 2 * 8] = gb2_xi * b2_eta;
        dN[base + 1 + 2 * 8] = b2_xi  * gb2_eta;
    }
}

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-2] Read coords — n×2 double matrix (column-major)
        // [FIX-F] np is std::size_t
        // -----------------------------------------------------------------------
        const TypedArray<double> coordArr = inputs[0];
        const auto dims                   = coordArr.getDimensions();
        const std::size_t np              = dims[0];

        // [FIX-B] Copy into contiguous vector for raw pointer access in kernel
        std::vector<double> coordsVec(coordArr.begin(), coordArr.end());

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-4] Allocate flat output buffer: 2 × 9 × np elements
        //                Fill via kernel then copy into TypedArray
        // -----------------------------------------------------------------------
        std::vector<double> dN_vec(2 * 9 * np, 0.0);

        computeDerivatives(dN_vec.data(), coordsVec.data(), np);

        // -----------------------------------------------------------------------
        // [NEW-1] Build 3D TypedArray output: {2, 9, np}
        //         Column-major flat layout matches dN[row + 2*col + 18*i] exactly
        // -----------------------------------------------------------------------
        TypedArray<double> dN_out =
            factory.createArray<double>({2, 9, np});
        std::copy(dN_vec.begin(), dN_vec.end(), dN_out.begin());

        outputs[0] = std::move(dN_out);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("Quad9Deriv:inputError",
                       "One input required: list (n x 2 array).");

        if (outputs.size() > 1)
            throwError("Quad9Deriv:outputError",
                       "One output required: dN (2 x 9 x n array).");

        // Must be a real double array
        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("Quad9Deriv:inputError",
                       "Input must be a real double n x 2 array.");

        // Must have exactly 2 columns [NEW-2]
        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 2)
            throwError("Quad9Deriv:inputError",
                       "Input must be a real double n x 2 array.");
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
