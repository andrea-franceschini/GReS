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
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
//
// [FIX-B] TypedArray<T>::operator[] returns a proxy — cannot take address.
//         mxToPolygon() now takes a TypedArray<double> by const-ref and
//         reads dimensions/elements via the C++ MEX API instead of
//         mxGetM / mxGetPr.  polygonToMx() replaced by polygonToArray()
//         returning a TypedArray<double> built with ArrayFactory.
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         String arguments use factory.createCharArray().
//         Note: bool IS arithmetic so factory.createScalar<bool>(valid)
//         is legal and replaces mxCreateLogicalScalar().
//
// [FIX-F] mwSize is declared in mex.h — NOT included by mex.hpp.
//         Every mwSize replaced with std::size_t (#include <cstddef>).
//
// [NEW-1] mxToPolygon() and polygonToMx() refactored as private helpers
//         of MexFunction so they have access to ArrayFactory and the
//         C++ MEX API, removing all legacy mxGetM / mxGetPr / mxGetN /
//         mxCreateDoubleMatrix calls.
//
// [NEW-2] mxCreateLogicalScalar(valid) replaced with
//         factory.createScalar<bool>(valid).
//
// [NEW-3] Empty polygon output (invalid clip) written as a 0 x 2 array
//         via factory.createArray<double>({0, 2}) — matches the original
//         polygonToMx(Polygon{}) which produced a 0x2 matrix.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "polygonClip.hpp"

#include <cstddef>   // std::size_t  [FIX-F]
#include <vector>
#include <string>

using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    //------------------------------------------------------------------------------------
    // [NEW-1][FIX-B][FIX-F] Replace mxToPolygon(const mxArray*)
    //   TypedArray<double> gives dimensions and element access without mxGetPr/mxGetM.
    //------------------------------------------------------------------------------------
    Polygon arrayToPolygon(const TypedArray<double>& arr)
    {
        const auto dims    = arr.getDimensions();
        const std::size_t N = dims[0];           // [FIX-F]

        Polygon poly(N);
        for (std::size_t i = 0; i < N; ++i)
            poly[i] = { arr[i], arr[i + N] };    // col-major: col0=x, col1=y

        return poly;
    }

    //------------------------------------------------------------------------------------
    // [NEW-1][FIX-B] Replace polygonToMx(const Polygon&) -> mxArray*
    //   Returns a TypedArray<double> (N x 2) built with ArrayFactory.
    // [NEW-3] Handles empty Polygon (N==0) → 0 x 2 array.
    //------------------------------------------------------------------------------------
    TypedArray<double> polygonToArray(const Polygon& poly)
    {
        const std::size_t N = poly.size();       // [FIX-F]

        TypedArray<double> arr = factory.createArray<double>({N, 2});

        for (std::size_t i = 0; i < N; ++i) {
            arr[i]     = poly[i][0];             // col 0: x
            arr[i + N] = poly[i][1];             // col 1: y
        }

        return arr;
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // [FIX-B] Convert inputs to Polygon via C++ API helper
        // -----------------------------------------------------------------------
        const TypedArray<double> polyArr    = inputs[0];
        const TypedArray<double> clipArr    = inputs[1];

        Polygon poly    = arrayToPolygon(polyArr);
        Polygon clipper = arrayToPolygon(clipArr);

        // -----------------------------------------------------------------------
        // Clip — logic unchanged
        // -----------------------------------------------------------------------
        bool valid = false;
        if (poly.size() >= 3 && clipper.size() >= 3)
            valid = clipPolygon(poly, clipper);

        // -----------------------------------------------------------------------
        // [NEW-2] Pack outputs
        //   polyOut : clipped polygon (or 0x2 if invalid)  [NEW-3]
        //   isValid : bool scalar replacing mxCreateLogicalScalar
        // -----------------------------------------------------------------------
        outputs[0] = polygonToArray(valid ? poly : Polygon{});
        outputs[1] = factory.createScalar<bool>(valid);   // [NEW-2]
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 2)
            throwError("polygonClip:inputError",
                       "Usage: [polyOut, isValid] = polygonClip(poly, clipPoly).");

        if (outputs.size() > 2)
            throwError("polygonClip:outputError",
                       "Two outputs required.");

        // Both inputs must be real double N x 2 matrices
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
