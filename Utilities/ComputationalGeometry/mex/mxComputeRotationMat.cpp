//----------------------------------------------------------------------------------------
// computeRot_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   R = computeRot(n)
//     n : 3-element real double vector (normal)
//     R : 3 x 3 real double rotation matrix (columns: n, m1, m2)
//
// Build command:
//   mex -R2018a computeRot_wrap.cpp
//
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() in the gateway for all validation
//         and post-computation checks.
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         n input copied into std::vector<double>; elements read into
//         std::array<double,3> for the math code.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         All string arguments use factory.createCharArray() instead.
//
// [NEW-1] normalize() called mexErrMsgTxt() directly from a pure math
//         function — same anti-pattern as inv3x3() in the Jacobian file.
//         Fixed by:
//           - Removing mexErrMsgTxt from normalize()
//           - Having normalize() throw std::runtime_error instead
//           - Catching it in operator()() and re-throwing via throwError()
//         This keeps math helpers MATLAB-free and ensures RAII unwinds.
//
// [NEW-2] mxGetPr() / mxCreateDoubleMatrix(3,3) replaced with
//         TypedArray<double> copy (input) and
//         factory.createArray<double>({3, 3}) (output).
//
// [NEW-3] Output filled via flat column-major index i + j*3, identical
//         to the original R_ptr[i + j*3] = Rcols[j][i] loop.
//
// [NEW-4] All pure math helpers (cross, dot, norm, normalize, det3x3)
//         are kept unchanged in logic — only mexErrMsgTxt removed from
//         normalize() and replaced with std::runtime_error.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <array>
#include <cmath>
#include <stdexcept>   // std::runtime_error  [NEW-1]
#include <vector>
#include <string>

using namespace matlab::data;
using matlab::mex::ArgumentList;

//========================================================================================
// Pure math helpers — NO MATLAB API dependency
// [NEW-1] normalize() throws std::runtime_error instead of mexErrMsgTxt
//========================================================================================

inline std::array<double,3> cross(const std::array<double,3>& a,
                                  const std::array<double,3>& b)
{
    return { a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] };
}

inline double dot(const std::array<double,3>& a, const std::array<double,3>& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double norm(const std::array<double,3>& a)
{
    return std::sqrt(dot(a, a));
}

// [NEW-1] throws std::runtime_error — keeps function MATLAB-free,
//         lets RAII destructors unwind on error
inline std::array<double,3> normalize(const std::array<double,3>& a)
{
    const double n = norm(a);
    if (n < 1e-16)
        throw std::runtime_error("Zero-length vector cannot be normalized.");
    return { a[0]/n, a[1]/n, a[2]/n };
}

inline double det3x3(const std::array<std::array<double,3>,3>& R)
{
    const auto& a = R[0]; const auto& b = R[1]; const auto& c = R[2];
    return a[0]*(b[1]*c[2] - b[2]*c[1])
         - a[1]*(b[0]*c[2] - b[2]*c[0])
         + a[2]*(b[0]*c[1] - b[1]*c[0]);
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
        // [FIX-B][NEW-2] Read n — 3-element double vector
        // -----------------------------------------------------------------------
        const TypedArray<double> nArr = inputs[0];

        // Copy into std::array for the math code
        std::array<double,3> n = { nArr[0], nArr[1], nArr[2] };

        // -----------------------------------------------------------------------
        // Compute rotation matrix
        // [NEW-1] normalize() may throw std::runtime_error — caught here and
        //         re-thrown as MATLAB exception so RAII destructors run cleanly
        // -----------------------------------------------------------------------
        try {
            n = normalize(n);

            // Pick a vector not parallel to n
            std::array<double,3> tmp;
            if (std::fabs(n[0]) < 0.9) tmp = {1.0, 0.0, 0.0};
            else                        tmp = {0.0, 1.0, 0.0};

            // First tangent: orthogonalize tmp against n
            const double ndot = dot(tmp, n);
            std::array<double,3> m1 = { tmp[0] - ndot*n[0],
                                        tmp[1] - ndot*n[1],
                                        tmp[2] - ndot*n[2] };
            m1 = normalize(m1);

            // Second tangent: orthogonal to both
            std::array<double,3> m2 = cross(n, m1);
            m2 = normalize(m2);

            // Assemble rotation matrix (columns: n, m1, m2)
            std::array<std::array<double,3>,3> Rcols = { n, m1, m2 };

            // Enforce right-handed system (det = +1)
            double detR = det3x3(Rcols);
            if (detR < 0.0) {
                m1[0] = -m1[0]; m1[1] = -m1[1]; m1[2] = -m1[2];
                Rcols = { n, m1, m2 };
                detR  = det3x3(Rcols);
            }

            if (std::fabs(detR - 1.0) > 1e-12)
                throwError("computeRot:notOrthogonal",
                           "Rotation matrix is not orthogonal to machine precision.");

            // -------------------------------------------------------------------
            // [NEW-2][NEW-3] Build 3x3 output array — column-major
            //   R_out[i + j*3] = Rcols[j][i]   (same as original)
            // -------------------------------------------------------------------
            TypedArray<double> R_out = factory.createArray<double>({3, 3});

            std::size_t k = 0;
            for (int j = 0; j < 3; ++j)
                for (int i = 0; i < 3; ++i, ++k)
                    R_out[k] = Rcols[j][i];

            outputs[0] = std::move(R_out);

        } catch (const std::runtime_error& e) {
            // [NEW-1] Re-throw as MATLAB exception — RAII unwinds cleanly
            throwError("computeRot:mathError", e.what());
        }
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("computeRot:inputError",
                       "Usage: R = computeRot(n).");

        if (outputs.size() > 1)
            throwError("computeRot:outputError",
                       "One output required.");

        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("computeRot:inputError",
                       "Input n must be a 3-element real vector.");

        if (inputs[0].getNumberOfElements() != 3)
            throwError("computeRot:inputError",
                       "Input n must be a 3-element real vector.");
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
