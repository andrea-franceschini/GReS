//----------------------------------------------------------------------------------------
// mxPolygonClip.cpp
// MEX gateway — MathWorks C++ MEX API (R2018a+)
//
// MATLAB signature:
//   [polyOut, isValid] = polygonClip(poly, clipPoly)
//     poly     : N x 2 real double
//     clipPoly : M x 2 real double
//     polyOut  : K x 2 real double (clipped polygon, or 0x2 if invalid)
//     isValid  : logical scalar
//
// Build command:
//   mex -R2018a mxPolygonClip.cpp
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "src/polygonClip.hpp"
#include "src/MexHelper.hpp"

#include <cstddef>
#include <vector>
#include <algorithm>

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace mexhelper;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    void throwErr(const std::string& id, const std::string& msg) {
        throwError(this, factory, id, msg);
    }

    Polygon arrayToPolygon(const TypedArray<double>& arr)
    {
        const auto dims     = arr.getDimensions();
        const std::size_t N = dims[0];
        std::vector<double> buf(arr.begin(), arr.end());
        Polygon poly(N);
        for (std::size_t i = 0; i < N; ++i)
            poly[i] = { buf[i], buf[i + N] };
        return poly;
    }

    TypedArray<double> polygonToArray(const Polygon& poly)
    {
        const std::size_t N = poly.size();
        std::vector<double> buf(N * 2, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            buf[i]     = poly[i][0];
            buf[i + N] = poly[i][1];
        }
        TypedArray<double> arr = factory.createArray<double>({N, 2});
        std::copy(buf.begin(), buf.end(), arr.begin());
        return arr;
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        if (inputs.size() != 2)
            throwErr("polygonClip:inputError",
                     "Usage: [polyOut, isValid] = polygonClip(poly, clipPoly).");

        if (outputs.size() > 2)
            throwErr("polygonClip:outputError", "Two outputs required.");

        for (std::size_t i = 0; i < 2; ++i) {
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwErr("polygonClip:inputError",
                         "Inputs must be N x 2 real double matrices.");
            const auto dims = inputs[i].getDimensions();
            if (dims.size() < 2 || dims[1] != 2)
                throwErr("polygonClip:inputError",
                         "Inputs must be N x 2 real double matrices.");
        }

        Polygon poly    = arrayToPolygon(inputs[0]);
        Polygon clipper = arrayToPolygon(inputs[1]);

        bool valid = false;
        if (poly.size() >= 3 && clipper.size() >= 3)
            valid = clipPolygon(poly, clipper);

        outputs[0] = valid ? polygonToArray(poly)
                           : factory.createArray<double>({0, 2});
        outputs[1] = factory.createScalar<bool>(valid);
    }
};

//----------------------------------------------------------------------------------------