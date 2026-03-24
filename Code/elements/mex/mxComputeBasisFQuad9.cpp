//----------------------------------------------------------------------------------------
// computeQuad9Basis_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   N = computeQuad9Basis(coords)   % coords: n×2, N: n×9
//
// Build command:
//   mex -R2018a computeQuad9Basis_wrap.cpp
//
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-F] mwSize is declared in mex.h — which is NOT included by the
//         pure C++ MEX API (mex.hpp). Using it caused "unknown type name
//         'mwSize'" errors at every use site.
//         Fixed by replacing every mwSize with std::size_t throughout,
//         both in the computational kernel and in the gateway.
//         std::size_t is the correct portable equivalent; it is the same
//         underlying type as mwSize on all 64-bit platforms.
//
// All previously documented fixes (FIX-A/B/C/E, NEW-1..4) are retained.
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
// Pure computational kernel — no MATLAB API dependency.
// [FIX-F] mwSize replaced with std::size_t
//
// coords layout: column-major, np rows × 2 cols
//   xi  = coords[i]       (col 0)
//   eta = coords[i + np]  (col 1)
// N layout: column-major, np rows × 9 cols
//   N[i + np*j] = shape function j at point i
//----------------------------------------------------------------------------------------
static void computeBasis(double* N, const double* coords, std::size_t np)
{
    for (std::size_t i = 0; i < np; ++i) {
        const double xi  = coords[i];         // col 0
        const double eta = coords[i + np];    // col 1

        // 1D quadratic basis functions
        const double b1_xi  = 0.5 * xi  * (xi  - 1.0);
        const double b2_xi  = 1.0 - xi  * xi;
        const double b3_xi  = 0.5 * xi  * (xi  + 1.0);

        const double b1_eta = 0.5 * eta * (eta - 1.0);
        const double b2_eta = 1.0 - eta * eta;
        const double b3_eta = 0.5 * eta * (eta + 1.0);

        // Tensor-product shape functions for Quad9
        N[i + np * 0] = b1_xi * b1_eta;  // Node 1
        N[i + np * 1] = b3_xi * b1_eta;  // Node 2
        N[i + np * 2] = b3_xi * b3_eta;  // Node 3
        N[i + np * 3] = b1_xi * b3_eta;  // Node 4
        N[i + np * 4] = b2_xi * b1_eta;  // Node 5
        N[i + np * 5] = b3_xi * b2_eta;  // Node 6
        N[i + np * 6] = b2_xi * b3_eta;  // Node 7
        N[i + np * 7] = b1_xi * b2_eta;  // Node 8
        N[i + np * 8] = b2_xi * b2_eta;  // Node 9
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
        // Read coords — n×2 double matrix (column-major)
        // [FIX-F] mwSize → std::size_t
        // -----------------------------------------------------------------------
        const TypedArray<double> coordArr = inputs[0];
        const auto dims                   = coordArr.getDimensions();
        const std::size_t np              = dims[0];   // [FIX-F]

        // Copy into contiguous vector for raw pointer access in kernel
        std::vector<double> coordsVec(coordArr.begin(), coordArr.end());

        // -----------------------------------------------------------------------
        // Allocate output work buffer, fill via kernel, copy into TypedArray
        // -----------------------------------------------------------------------
        std::vector<double> N_vec(np * 9, 0.0);

        computeBasis(N_vec.data(), coordsVec.data(), np);

        TypedArray<double> N_out = factory.createArray<double>({np, 9});
        std::copy(N_vec.begin(), N_vec.end(), N_out.begin());

        outputs[0] = std::move(N_out);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("Quad9:inputError",
                       "One input required: coordList (n x 2 array).");

        if (outputs.size() > 1)
            throwError("Quad9:outputError",
                       "One output required: N (n x 9 array).");

        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("Quad9:inputError",
                       "Input must be a real double n x 2 array.");

        const auto dims = inputs[0].getDimensions();
        if (dims.size() < 2 || dims[1] != 2)
            throwError("Quad9:inputError",
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
