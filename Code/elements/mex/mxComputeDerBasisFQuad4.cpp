//----------------------------------------------------------------------------------------
// computeQuad4Derivatives_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   dN = computeQuad4Derivatives(coords)   % coords: n×2, dN: 2×4×n
//
// Build command:
//   mex -R2018a computeQuad4Derivatives_wrap.cpp
//
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
//         Generic IDs assigned since original used no-ID mexErrMsgTxt.
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
//         Affected sites:
//           - computeQuad4Derivatives() parameter types
//           - both loop variables (i, j)
//           - np declaration in operator()()
//
// [NEW-1] mxGetPr() / mxCreateNumericArray(3, {2,4,np}) replaced with
//         factory.createArray<double>({2, 4, np}) and std::copy.
//
// [NEW-2] mxGetM() / mxGetN() replaced with getDimensions()[0] / [1].
//
// [NEW-3] The pure computational function computeQuad4Derivatives() is
//         kept unchanged in logic — no MATLAB API dependency in kernel.
//
// [NEW-4] 3D output indexing: dN[row + 2*col + 8*i]
//         factory.createArray<double>({2, 4, np}) produces a column-major
//         3D array with stride 2*4 = 8 between point slices — identical
//         to the original mxCreateNumericArray layout.
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
// [FIX-F] mwSize → std::size_t throughout
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
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-2] Read coords — n×2 double matrix (column-major)
        // [FIX-F] np is std::size_t, no cast needed
        // -----------------------------------------------------------------------
        const TypedArray<double> coordArr = inputs[0];
        const auto dims                   = coordArr.getDimensions();
        const std::size_t np              = dims[0];   // [FIX-F]

        // [FIX-B] Copy into contiguous vector for raw pointer access in kernel
        std::vector<double> coordsVec(coordArr.begin(), coordArr.end());

        // -----------------------------------------------------------------------
        // [FIX-B][NEW-4] Allocate flat output buffer: 2 × 4 × np elements
        //                Fill via kernel, then copy into TypedArray
        // -----------------------------------------------------------------------
        std::vector<double> dN_vec(2 * 4 * np, 0.0);

        computeQuad4Derivatives(dN_vec.data(), coordsVec.data(), np);

        // -----------------------------------------------------------------------
        // [NEW-1] Build 3D TypedArray output: {2, 4, np}
        //         Column-major flat layout matches dN[row + 2*col + 8*i] exactly
        // -----------------------------------------------------------------------
        TypedArray<double> dN_out =
            factory.createArray<double>({2, 4, np});
        std::copy(dN_vec.begin(), dN_vec.end(), dN_out.begin());

        outputs[0] = std::move(dN_out);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
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
