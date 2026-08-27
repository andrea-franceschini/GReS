//----------------------------------------------------------------------------------------
// computeQuad9Derivatives_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   dN = computeQuad9Derivatives(coords)   % coords: n×2, dN: 2×9×n
//
//----------------------------------------------------------------------------------------

#include <cstdint>
#include "mex.hpp"
#include "mexAdapter.hpp"

#include <cstddef>   // std::size_t 
#include <vector>
#include <string>
#include <algorithm> // std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
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
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        const TypedArray<double> coordArr = inputs[0];
        const auto dims                   = coordArr.getDimensions();
        const std::size_t np              = dims[0];

        // Copy into contiguous vector for raw pointer access in kernel
        std::vector<double> coordsVec(coordArr.begin(), coordArr.end());

        // -----------------------------------------------------------------------
        // Fill via kernel then copy into TypedArray
        // -----------------------------------------------------------------------
        std::vector<double> dN_vec(2 * 9 * np, 0.0);

        computeDerivatives(dN_vec.data(), coordsVec.data(), np);

        // -----------------------------------------------------------------------
        // Column-major flat layout matches dN[row + 2*col + 18*i] exactly
        // -----------------------------------------------------------------------
        TypedArray<double> dN_out =
            factory.createArray<double>({2, 9, np});
        std::copy(dN_vec.begin(), dN_vec.end(), dN_out.begin());

        outputs[0] = std::move(dN_out);
    }

private:

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
