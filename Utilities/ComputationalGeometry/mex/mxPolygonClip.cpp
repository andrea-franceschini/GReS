//----------------------------------------------------------------------------------------
// polygonClip_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [polyOut, isValid] = polygonClip(poly, clipPoly)
//     poly     : N x 2 real double
//     clipPoly : M x 2 real double
//     polyOut  : K x 2 real double (clipped polygon, or 0x2 if invalid)
//     isValid  : logical scalar
//
// Build command:
//   mex -R2018a polygonClip_wrap.cpp
//
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-G] TypedArray flat-index subscripting causes runtime error
//         "Not enough indices provided" on 2D arrays.
//         TypedArray::operator[](size_t) on a {N,2} array does NOT do
//         flat/linear indexing — it expects two subscripts.
//         TypedArray iterators DO traverse in flat column-major order.
//
//         arrayToPolygon(): instead of arr[i] / arr[i+N] on a 2D
//         TypedArray, copy the whole array into a std::vector<double>
//         first and index into that.
//
//         polygonToArray(): instead of arr[i] = x / arr[i+N] = y on a
//         2D TypedArray, fill a std::vector<double> column-major and
//         std::copy into the TypedArray via iterators.
//
// All previously documented fixes (FIX-A/B/C/E/F, NEW-1..3) are retained.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "polygonClip.hpp"

#include <cstddef>   // std::size_t
#include <vector>
#include <string>
#include <algorithm> // std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    //------------------------------------------------------------------------------------
    // [FIX-G] arrayToPolygon
    //   Copy TypedArray into std::vector first, then index into the vector.
    //   Avoids flat-index subscripting on the 2D TypedArray.
    //------------------------------------------------------------------------------------
    Polygon arrayToPolygon(const TypedArray<double>& arr)
    {
        const auto dims     = arr.getDimensions();
        const std::size_t N = dims[0];

        // Copy into flat vector — iterators traverse column-major order
        std::vector<double> buf(arr.begin(), arr.end());

        Polygon poly(N);
        for (std::size_t i = 0; i < N; ++i)
            poly[i] = { buf[i], buf[i + N] };   // col 0 = x, col 1 = y

        return poly;
    }

    //------------------------------------------------------------------------------------
    // [FIX-G] polygonToArray
    //   Fill a flat std::vector column-major, then std::copy into TypedArray.
    //   Avoids flat-index subscripting on the 2D TypedArray.
    //------------------------------------------------------------------------------------
    TypedArray<double> polygonToArray(const Polygon& poly)
    {
        const std::size_t N = poly.size();

        // Fill flat column-major buffer: col 0 = x, col 1 = y
        std::vector<double> buf(N * 2, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            buf[i]     = poly[i][0];   // col 0: x
            buf[i + N] = poly[i][1];   // col 1: y
        }

        TypedArray<double> arr = factory.createArray<double>({N, 2});
        std::copy(buf.begin(), buf.end(), arr.begin());
        return arr;
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        Polygon poly    = arrayToPolygon(inputs[0]);
        Polygon clipper = arrayToPolygon(inputs[1]);

        bool valid = false;
        if (poly.size() >= 3 && clipper.size() >= 3)
            valid = clipPolygon(poly, clipper);

        outputs[0] = polygonToArray(valid ? poly : Polygon{});
        outputs[1] = factory.createScalar<bool>(valid);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 2)
            throwError("polygonClip:inputError",
                       "Usage: [polyOut, isValid] = polygonClip(poly, clipPoly).");

        if (outputs.size() > 2)
            throwError("polygonClip:outputError",
                       "Two outputs required.");

        for (std::size_t i = 0; i < 2; ++i) {
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("polygonClip:inputError",
                           "Inputs must be N x 2 real double matrices.");

            const auto dims = inputs[i].getDimensions();
            if (dims.size() < 2 || dims[1] != 2)
                throwError("polygonClip:inputError",
                           "Inputs must be N x 2 real double matrices.");
        }
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
