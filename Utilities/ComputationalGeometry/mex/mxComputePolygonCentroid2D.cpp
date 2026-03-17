//----------------------------------------------------------------------------------------
// computePolygonCentroid2D_CCW_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   centroid = computePolygonCentroid2D_CCW(points)
//     points   : N x 2 real double (CCW ordered)
//     centroid : 1 x 2 real double [Cx, Cy]
//
// Build command:
//   mex -R2018a computePolygonCentroid2D_CCW_wrap.cpp
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
//         pointer arithmetic in the centroid computation loop.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         All string arguments use factory.createCharArray() instead.
//
// [FIX-F] mwSize is declared in mex.h — NOT included by mex.hpp.
//         Every mwSize replaced with std::size_t (#include <cstddef>).
//         Affected sites: N, dim, loop variable i, index j.
//
// [NEW-1] mxGetM() / mxGetN() replaced with getDimensions()[0] / [1].
//
// [NEW-2] mxCreateDoubleMatrix(1,2) + mxGetPr replaced with
//         factory.createArray<double>({1, 2}) and element assignment.
//
// [NEW-3] Zero-area check throws via throwError() — the std::vector<double>
//         P is live at that point so the old mexErrMsgIdAndTxt longjmp
//         would have leaked it; a C++ exception unwinds it cleanly.
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
        // [FIX-B][NEW-1] Read points — N x 2 double matrix (column-major)
        // [FIX-F] N, dim are std::size_t
        // -----------------------------------------------------------------------
        const TypedArray<double> Parr = inputs[0];
        const auto dims               = Parr.getDimensions();   // [NEW-1]
        const std::size_t N           = dims[0];                // [FIX-F]

        // [FIX-B] Copy into contiguous vector for raw pointer arithmetic
        std::vector<double> P(Parr.begin(), Parr.end());

        // -----------------------------------------------------------------------
        // Centroid computation (shoelace / Green's theorem)
        // [FIX-F] loop variables i, j are std::size_t
        // -----------------------------------------------------------------------
        double Cx = 0.0;
        double Cy = 0.0;
        double A  = 0.0;   // signed area × 2

        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t j = (i + 1) % N;   // [FIX-F]

            const double xi = P[i];
            const double yi = P[i + N];

            const double xj = P[j];
            const double yj = P[j + N];

            const double cross = xi*yj - xj*yi;

            A  += cross;
            Cx += (xi + xj) * cross;
            Cy += (yi + yj) * cross;
        }

        // [NEW-3] Zero-area check via throwError() — C++ exception unwinds P
        if (A == 0.0)
            throwError("Centroid2D:zeroArea", "Zero-area polygon.");

        A  *= 0.5;
        Cx /= (6.0 * A);
        Cy /= (6.0 * A);

        // -----------------------------------------------------------------------
        // [NEW-2] Build output: 1 x 2 double array
        //         factory.createArray replaces mxCreateDoubleMatrix + mxGetPr
        // -----------------------------------------------------------------------
        TypedArray<double> out = factory.createArray<double>({1, 2});
        out[0] = Cx;
        out[1] = Cy;

        outputs[0] = std::move(out);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("Centroid2D:inputError",
                       "Usage: centroid = computePolygonCentroid2D_CCW(points).");

        if (outputs.size() > 1)
            throwError("Centroid2D:outputError",
                       "One output required.");

        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("Centroid2D:inputError",
                       "points must be a real double array.");

        // Must be N x 2 [NEW-1]
        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 2)
            throwError("Centroid2D:inputError",
                       "points must be N x 2.");

        // At least 3 points required
        if (dims[0] < 3)
            throwError("Centroid2D:inputError",
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
