#pragma once

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "PolygonGeometry.hpp"

#include <cstddef>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace matlab::data;
using matlab::mex::ArgumentList;

namespace mexhelper {

inline void throwError(matlab::mex::Function* fn,
                       ArrayFactory& factory,
                       const std::string& id,
                       const std::string& msg)
{
    fn->getEngine()->feval(u"error", 0,
        std::vector<Array>{
            factory.createCharArray(id),
            factory.createCharArray(msg)
        });
}

inline void requireArgs(matlab::mex::Function* fn,
                        ArrayFactory& factory,
                        bool cond,
                        const std::string& id,
                        const std::string& msg)
{
    if (!cond) throwError(fn, factory, id, msg);
}

inline std::vector<double> typedArrayToVector(const TypedArray<double>& arr)
{
    return std::vector<double>(arr.begin(), arr.end());
}

inline TypedArray<double> vectorToArray(ArrayFactory& factory,
                                        const std::vector<double>& v,
                                        std::size_t rows,
                                        std::size_t cols)
{
    TypedArray<double> arr = factory.createArray<double>({rows, cols});
    std::copy(v.begin(), v.end(), arr.begin());
    return arr;
}

inline void validateDoubleMatrix(matlab::mex::Function* fn,
                                 ArrayFactory& factory,
                                 const Array& arr,
                                 const std::string& id,
                                 const std::string& ctx)
{
    requireArgs(fn, factory, arr.getType() == ArrayType::DOUBLE, id,
                ctx + " must be a real double matrix.");
    // Fix: isComplex() belongs to TypedArray, but we can check the Array directly or via property
    // For R2018a+ API, complex arrays have complex types (e.g. ArrayType::COMPLEX_DOUBLE)
    // The previous ArrayType::DOUBLE check inherently guarantees it is strictly real double.
}

inline std::size_t getRows(const Array& arr) { return arr.getDimensions()[0]; }
inline std::size_t getCols(const Array& arr) { return arr.getDimensions()[1]; }
inline std::size_t numel(const Array& arr) {
    std::size_t n = 1;
    for (auto d : arr.getDimensions()) n *= d;
    return n;
}

inline polygeom::BatchInput parseBatchInput(matlab::mex::Function* fn,
                                            ArrayFactory& factory,
                                            const Array& Pflat_arr,
                                            const Array& nVert_arr,
                                            std::vector<double>& Pbuf,
                                            std::vector<double>& nVertBuf)
{
    validateDoubleMatrix(fn, factory, Pflat_arr, "PolygonGeometry:input", "Pflat");
    validateDoubleMatrix(fn, factory, nVert_arr, "PolygonGeometry:input", "nVert");

    std::size_t nPts = getRows(Pflat_arr);
    std::size_t dim  = getCols(Pflat_arr);
    requireArgs(fn, factory, dim == 2 || dim == 3,
                "PolygonGeometry:input", "Pflat must be Ntot x 2 or Ntot x 3.");

    std::size_t nPoly = numel(nVert_arr);
    requireArgs(fn, factory, nPoly >= 1,
                "PolygonGeometry:input", "nVert must contain at least one polygon.");

    Pbuf     = typedArrayToVector(Pflat_arr);
    nVertBuf = typedArrayToVector(nVert_arr);

    double sum = 0.0;
    for (std::size_t i = 0; i < nPoly; ++i) {
        requireArgs(fn, factory,
                    polygeom::isIntegerValued(nVertBuf[i]) && nVertBuf[i] >= 3.0,
                    "PolygonGeometry:input", "Each nVert entry must be an integer >= 3.");
        sum += nVertBuf[i];
    }
    requireArgs(fn, factory,
                static_cast<std::size_t>(std::llround(sum)) == nPts,
                "PolygonGeometry:input", "sum(nVert) must equal size(Pflat,1).");

    polygeom::BatchInput in;
    in.P     = Pbuf.data();
    in.nPts  = nPts;
    in.dim   = dim;
    in.nVert = nVertBuf.data();
    in.nPoly = nPoly;
    return in;
}

// Fix: Remove 'const' from ArgumentList as the MATLAB implementation doesn't support it for size/operator[]
inline bool isLikelyBatchCall(ArgumentList& inputs)
{
    if (inputs.size() < 2) return false;
    if (inputs[1].getType() == ArrayType::CHAR ||
        inputs[1].getType() == ArrayType::MATLAB_STRING) return false;
    bool secondIsVector = (getRows(inputs[1]) == 1 || getCols(inputs[1]) == 1);
    bool firstLooksMatrix = (getCols(inputs[0]) == 2 || getCols(inputs[0]) == 3);
    bool secondHasMany = (numel(inputs[1]) > 1);
    return firstLooksMatrix && secondIsVector && secondHasMany;
}

inline bool isCharArray(const Array& arr) {
    return arr.getType() == ArrayType::CHAR ||
           arr.getType() == ArrayType::MATLAB_STRING;
}

inline std::string charArrayToString(const Array& arr) {
    if (arr.getType() == ArrayType::MATLAB_STRING)
        return std::string(CharArray(arr).toAscii());
    return std::string(CharArray(arr).toAscii());
}

} // namespace mexhelper