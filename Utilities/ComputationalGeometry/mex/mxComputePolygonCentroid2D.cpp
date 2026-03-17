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
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-G] TypedArray flat-index subscripting causes runtime error
//         "Not enough indices provided" on 2D arrays.
//
//         out {1, 2}:
//           out[0] and out[1] were single flat indices on a 2D TypedArray.
//           Fixed: write into a std::vector<double>{Cx, Cy} and
//           std::copy into the TypedArray via iterators.
//
// All previously documented fixes (FIX-A/B/C/E/F, NEW-1..3) are retained.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cstddef>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm> // std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        const TypedArray<double> Parr = inputs[0];
        const auto dims               = Parr.getDimensions();
        const std::size_t N           = dims[0];

        std::vector<double> P(Parr.begin(), Parr.end());

        // Shoelace / Green's theorem
        double Cx = 0.0, Cy = 0.0, A = 0.0;

        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t j = (i + 1) % N;

            const double xi = P[i],     yi = P[i + N];
            const double xj = P[j],     yj = P[j + N];
            const double cr = xi*yj - xj*yi;

            A  += cr;
            Cx += (xi + xj) * cr;
            Cy += (yi + yj) * cr;
        }

        if (A == 0.0)
            throwError("Centroid2D:zeroArea", "Zero-area polygon.");

        A  *= 0.5;
        Cx /= (6.0 * A);
        Cy /= (6.0 * A);

        // -----------------------------------------------------------------------
        // [FIX-G] Write via std::vector + std::copy — not direct arr[0]/arr[1]
        // -----------------------------------------------------------------------
        std::vector<double> outBuf = { Cx, Cy };
        TypedArray<double> out = factory.createArray<double>({1, 2});
        std::copy(outBuf.begin(), outBuf.end(), out.begin());

        outputs[0] = std::move(out);
    }

private:

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

        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 2)
            throwError("Centroid2D:inputError",
                       "points must be N x 2.");

        if (dims[0] < 3)
            throwError("Centroid2D:inputError",
                       "At least 3 points required.");
    }

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
