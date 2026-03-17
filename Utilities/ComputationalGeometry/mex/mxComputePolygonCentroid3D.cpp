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
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-G] TypedArray flat-index subscripting causes runtime error
//         "Not enough indices provided" on 2D arrays.
//
//         out {1, 3}:
//           out[0], out[1], out[2] were single flat indices on a 2D
//           TypedArray. Fixed: write into std::vector<double>{C[0..2]}
//           and std::copy into the TypedArray via iterators.
//
// All previously documented fixes (FIX-A/B/C/E/F, NEW-1..4) are retained.
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

        const double P0[3] = { P[0], P[N], P[2*N] };

        double C[3]  = {0.0, 0.0, 0.0};
        double Atot  = 0.0;

        for (std::size_t i = 1; i < N - 1; ++i) {
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

        if (Atot == 0.0)
            throwError("Centroid3D:zeroArea", "Zero-area polygon.");

        C[0] /= Atot;
        C[1] /= Atot;
        C[2] /= Atot;

        // -----------------------------------------------------------------------
        // [FIX-G] Write via std::vector + std::copy — not direct arr[0..2]
        // -----------------------------------------------------------------------
        std::vector<double> outBuf = { C[0], C[1], C[2] };
        TypedArray<double> out = factory.createArray<double>({1, 3});
        std::copy(outBuf.begin(), outBuf.end(), out.begin());

        outputs[0] = std::move(out);
    }

private:

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

        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 3)
            throwError("Centroid3D:inputError",
                       "points must be N x 3.");

        if (dims[0] < 3)
            throwError("Centroid3D:inputError",
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
