//----------------------------------------------------------------------------------------
// computeQuad4Derivatives_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   dN = computeQuad4Derivatives(coords)   % coords: n×2, dN: 2×4×n
//
//----------------------------------------------------------------------------------------

#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cstddef>   // std::size_t  [FIX-F]
#include <vector>
#include <string>
#include <algorithm> // std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// Pure computational kernel — no MATLAB API dependency.
//
// coords layout: column-major, np rows × 2 cols
//   xi  = coords[i]       (col 0)
//   eta = coords[i + np]  (col 1)
//
// dN layout: column-major 3D, dims {2, 4, np}
//   dN[row + 2*col + 8*i]
//     row = 0 → d/dxi,  row = 1 → d/deta
//     col = 0..3 (node index)
//     i   = 0..np-1 (point index)
//----------------------------------------------------------------------------------------
static void computeQuad4Derivatives(double* dN, const double* coords, std::size_t np)
{
    for (std::size_t i = 0; i < np; ++i) {
        const double xi  = coords[i];         // col 0
        const double eta = coords[i + np];    // col 1

        const double dN_dxi[4] = {
            -0.25 * (1.0 - eta),   // Node 1
             0.25 * (1.0 - eta),   // Node 2
             0.25 * (1.0 + eta),   // Node 3
            -0.25 * (1.0 + eta)    // Node 4
        };

        const double dN_deta[4] = {
            -0.25 * (1.0 - xi),    // Node 1
            -0.25 * (1.0 + xi),    // Node 2
             0.25 * (1.0 + xi),    // Node 3
             0.25 * (1.0 - xi)     // Node 4
        };

        const std::size_t base = 8 * i;   // stride per point = 2 * 4

        for (std::size_t j = 0; j < 4; ++j) {
            dN[base + 0 + 2 * j] = dN_dxi [j];  // d/dxi
            dN[base + 1 + 2 * j] = dN_deta[j];  // d/deta
        }
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
        validateArguments(outputs, inputs);

        const TypedArray<double> coordArr = inputs[0];
        const auto dims                   = coordArr.getDimensions();
        const std::size_t np              = dims[0];   // [FIX-F]

        // Copy into contiguous vector for raw pointer access in kernel
        std::vector<double> coordsVec(coordArr.begin(), coordArr.end());

        // -----------------------------------------------------------------------
        // Fill via kernel, then copy into TypedArray
        // -----------------------------------------------------------------------
        std::vector<double> dN_vec(2 * 4 * np, 0.0);

        computeQuad4Derivatives(dN_vec.data(), coordsVec.data(), np);

        // -----------------------------------------------------------------------
        // Column-major flat layout matches dN[row + 2*col + 8*i] exactly
        // -----------------------------------------------------------------------
        TypedArray<double> dN_out =
            factory.createArray<double>({2, 4, np});
        std::copy(dN_vec.begin(), dN_vec.end(), dN_out.begin());

        outputs[0] = std::move(dN_out);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("Quad4Deriv:inputError",
                       "One input required: n x 2 coordinates.");

        if (outputs.size() > 1)
            throwError("Quad4Deriv:outputError",
                       "One output required: dN (2 x 4 x n array).");

        // Must be a real double array
        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("Quad4Deriv:inputError",
                       "Input must be a real double n x 2 array.");

        // Must have exactly 2 columns [NEW-2]
        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 2)
            throwError("Quad4Deriv:inputError",
                       "Input must be a real double n x 2 array.");
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
