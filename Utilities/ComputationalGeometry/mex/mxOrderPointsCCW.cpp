//----------------------------------------------------------------------------------------
// orderPointsCCW_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   idx = orderPointsCCW(points, normal)
//     points : N x 3 real double
//     normal : 3-element real double vector
//     idx    : N x 1 double (1-based sorted indices)
//
// Build command:
//   mex -R2018a orderPointsCCW_wrap.cpp
//
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgIdAndTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         Both input arrays (points, normal) copied into std::vector<double>;
//         .data() used for raw pointer arithmetic in the math code.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         All string arguments use factory.createCharArray() instead.
//
// [FIX-F] mwSize and mwIndex are declared in mex.h — NOT included by
//         mex.hpp. Both replaced with std::size_t throughout:
//           - AngleIndex::index field
//           - All loop variables (i)
//           - N, dim declarations
//           - Output index cast
//
// [NEW-1] mxGetDimensions() replaced with getDimensions() returning
//         std::vector<size_t>.
//
// [NEW-2] mxCreateDoubleMatrix + mxGetPr replaced with
//         factory.createArray<double>({N, 1}) and TypedArray iterator.
//
// [NEW-3] Degenerate centroid check kept — now throws via throwError()
//         instead of mexErrMsgIdAndTxt(), ensuring RAII unwinds cleanly.
//
// [NEW-4] Pure math helpers (dot3, cross3, norm3) are unchanged —
//         no MATLAB API dependency, operate on raw double* only.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cstddef>   // std::size_t  [FIX-F]
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// [FIX-F] mwIndex → std::size_t
//----------------------------------------------------------------------------------------
struct AngleIndex {
    double      angle;
    std::size_t index;   // [FIX-F]
};

//----------------------------------------------------------------------------------------
// [NEW-4] Pure math helpers — no MATLAB API dependency, unchanged logic
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
        // [FIX-B][NEW-1] Read input arrays into std::vector
        // [FIX-F] N and dim are std::size_t
        // -----------------------------------------------------------------------
        const TypedArray<double> pointsArr = inputs[0];
        const TypedArray<double> normalArr = inputs[1];

        const auto dims      = pointsArr.getDimensions();  // [NEW-1]
        const std::size_t N  = dims[0];                    // [FIX-F]
        const std::size_t dim = dims[1];                   // [FIX-F] (validated == 3)

        // [FIX-B] Copy into contiguous vectors for raw pointer arithmetic
        std::vector<double> points(pointsArr.begin(), pointsArr.end());
        std::vector<double> normal(normalArr.begin(), normalArr.end());

        // -----------------------------------------------------------------------
        // Compute centroid
        // [FIX-F] loop variable i is std::size_t
        // -----------------------------------------------------------------------
        double centroid[3] = {0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < N; ++i) {
            centroid[0] += points[i];
            centroid[1] += points[i +     N];
            centroid[2] += points[i + 2 * N];
        }
        centroid[0] /= static_cast<double>(N);
        centroid[1] /= static_cast<double>(N);
        centroid[2] /= static_cast<double>(N);

        // -----------------------------------------------------------------------
        // Reference vector v0 from point 0 to centroid
        // -----------------------------------------------------------------------
        double v0[3];
        v0[0] = centroid[0] - points[0];
        v0[1] = centroid[1] - points[N];
        v0[2] = centroid[2] - points[2 * N];

        const double nrm = norm3(v0);

        // [NEW-3] Degenerate check — throwError() instead of mexErrMsgIdAndTxt
        if (nrm == 0.0)
            throwError("orderPointsCCW:degenerate",
                       "First point coincides with centroid.");

        v0[0] /= nrm;
        v0[1] /= nrm;
        v0[2] /= nrm;

        // -----------------------------------------------------------------------
        // Compute angles relative to v0
        // [FIX-F] loop variable i is std::size_t
        // -----------------------------------------------------------------------
        std::vector<AngleIndex> angles(N);
        angles[0] = {0.0, 0};

        for (std::size_t i = 1; i < N; ++i) {
            double v[3];
            v[0] = centroid[0] - points[i];
            v[1] = centroid[1] - points[i +     N];
            v[2] = centroid[2] - points[i + 2 * N];

            const double dotv  = dot3(v, v0);

            double crossp[3];
            cross3(v, v0, crossp);

            const double detv = dot3(normal.data(), crossp);

            angles[i] = { std::atan2(detv, dotv), i };
        }

        // Sort CCW by angle
        std::sort(angles.begin(), angles.end(),
                  [](const AngleIndex& a, const AngleIndex& b) {
                      return a.angle < b.angle;
                  });

        // -----------------------------------------------------------------------
        // [NEW-2] Build output: N x 1 double, 1-based indices
        //         factory.createArray replaces mxCreateDoubleMatrix + mxGetPr
        // -----------------------------------------------------------------------
        TypedArray<double> out = factory.createArray<double>({N, 1});
        {
            auto it = out.begin();
            for (std::size_t i = 0; i < N; ++i, ++it)
                *it = static_cast<double>(angles[i].index + 1);  // 0-based → 1-based
        }

        outputs[0] = std::move(out);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 2)
            throwError("orderPointsCCW:nrhs",
                       "Two inputs required: points, normal.");

        if (outputs.size() > 1)
            throwError("orderPointsCCW:nlhs",
                       "One output required.");

        // points must be a real double array
        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("orderPointsCCW:points",
                       "points must be a real double array.");

        // normal must be a real double array
        if (inputs[1].getType() != ArrayType::DOUBLE)
            throwError("orderPointsCCW:normal",
                       "normal must be a real double array.");

        // points must be N x 3
        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 3)
            throwError("orderPointsCCW:points",
                       "points must be N x 3.");

        // normal must have exactly 3 elements
        if (inputs[1].getNumberOfElements() != 3)
            throwError("orderPointsCCW:normal",
                       "normal must have 3 elements.");
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
