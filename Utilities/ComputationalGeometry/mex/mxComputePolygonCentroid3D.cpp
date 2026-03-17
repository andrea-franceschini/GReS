//----------------------------------------------------------------------------------------
// computePolygonCentroidCCW_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   centroid = computePolygonCentroidCCW(points)
//     points   : N x 3 real double (CCW ordered, planar polygon in 3D)
//     centroid : 1 x 3 real double [Cx, Cy, Cz]
//
// Build command:
//   mex -R2018a computePolygonCentroidCCW_wrap.cpp
//
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
//         Generic IDs assigned since original used no-ID mexErrMsgTxt.
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         points copied into std::vector<double>; .data() used for raw
//         pointer arithmetic in the fan triangulation loop.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         All string arguments use factory.createCharArray() instead.
//
// [FIX-F] mwSize is declared in mex.h — NOT included by mex.hpp.
//         Every mwSize replaced with std::size_t (#include <cstddef>).
//         Affected sites: N, loop variable i.
//
// [NEW-1] mxGetM() / mxGetN() replaced with getDimensions()[0] / [1].
//
// [NEW-2] mxCreateDoubleMatrix(1,3) + mxGetPr replaced with
//         factory.createArray<double>({1, 3}) and element assignment.
//
// [NEW-3] Zero-area check throws via throwError() — std::vector<double> P
//         is live at that point so the old mexErrMsgTxt longjmp would have
//         leaked it; C++ exception unwinds it cleanly.
//
// [NEW-4] Pure math helpers (dot3, cross3, norm3) are unchanged —
//         no MATLAB API dependency.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cstddef>   // std::size_t  [FIX-F]
#include <vector>
#include <string>
#include <cmath>

using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// [NEW-4] Pure math helpers — no MATLAB API dependency, logic unchanged
//----------------------------------------------------------------------------------------
inline double dot3(const double* a, const double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void cross3(const double* a, const double* b, double* c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

inline double norm3(const double* a)
{
    return std::sqrt(dot3(a, a));
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
        // [FIX-B][NEW-1] Read points — N x 3 double matrix (column-major)
        // [FIX-F] N is std::size_t
        // -----------------------------------------------------------------------
        const TypedArray<double> Parr = inputs[0];
        const auto dims               = Parr.getDimensions();   // [NEW-1]
        const std::size_t N           = dims[0];                // [FIX-F]

        // [FIX-B] Copy into contiguous vector for raw pointer arithmetic
        std::vector<double> P(Parr.begin(), Parr.end());

        // -----------------------------------------------------------------------
        // Fan triangulation centroid (reference point = P[0])
        // [FIX-F] loop variable i is std::size_t
        // -----------------------------------------------------------------------
        const double P0[3] = { P[0], P[N], P[2*N] };

        double C[3]  = {0.0, 0.0, 0.0};
        double Atot  = 0.0;

        for (std::size_t i = 1; i < N - 1; ++i) {   // [FIX-F]
            const double v1[3] = {
                P[i]         - P0[0],
                P[i     + N] - P0[1],
                P[i + 2*N]   - P0[2]
            };

            const double v2[3] = {
                P[i+1]         - P0[0],
                P[i+1   + N]   - P0[1],
                P[i+1 + 2*N]   - P0[2]
            };

            double cp[3];
            cross3(v1, v2, cp);
            const double A = 0.5 * norm3(cp);

            const double Ct[3] = {
                (P0[0] + P[i]       + P[i+1])       / 3.0,
                (P0[1] + P[i + N]   + P[i+1 + N])   / 3.0,
                (P0[2] + P[i + 2*N] + P[i+1 + 2*N]) / 3.0
            };

            C[0] += A * Ct[0];
            C[1] += A * Ct[1];
            C[2] += A * Ct[2];
            Atot += A;
        }

        // [NEW-3] Zero-area check via throwError() — C++ exception unwinds P
        if (Atot == 0.0)
            throwError("Centroid3D:zeroArea", "Zero-area polygon.");

        C[0] /= Atot;
        C[1] /= Atot;
        C[2] /= Atot;

        // -----------------------------------------------------------------------
        // [NEW-2] Build output: 1 x 3 double array
        //         factory.createArray replaces mxCreateDoubleMatrix + mxGetPr
        // -----------------------------------------------------------------------
        TypedArray<double> out = factory.createArray<double>({1, 3});
        out[0] = C[0];
        out[1] = C[1];
        out[2] = C[2];

        outputs[0] = std::move(out);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("Centroid3D:inputError",
                       "Usage: centroid = computePolygonCentroidCCW(points).");

        if (outputs.size() > 1)
            throwError("Centroid3D:outputError",
                       "One output required.");

        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("Centroid3D:inputError",
                       "points must be a real double array.");

        // Must be N x 3  [NEW-1]
        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 3)
            throwError("Centroid3D:inputError",
                       "points must be N x 3.");

        // At least 3 points required
        if (dims[0] < 3)
            throwError("Centroid3D:inputError",
                       "At least 3 points required.");
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
