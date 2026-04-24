#include "mex.hpp"
#include "mexAdapter.hpp"
#include "src/MexHelper.hpp"

#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace mexhelper;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    enum class Mode { Geometry, Area, Centroid, Normal };

    void throwErr(const std::string& id, const std::string& msg) {
        throwError(this, factory, id, msg);
    }

    Mode parseMode(const Array& arr)
    {
        requireArgs(this, factory, isCharArray(arr),
                    "mxPolygonGeometry:input",
                    "Optional mode must be a character vector or string scalar.");
        std::string s = charArrayToString(arr);
        if (s == "geometry") return Mode::Geometry;
        if (s == "area")     return Mode::Area;
        if (s == "centroid") return Mode::Centroid;
        if (s == "normal")   return Mode::Normal;
        throwErr("mxPolygonGeometry:input",
                 "Mode must be one of: 'geometry', 'area', 'centroid', 'normal'.");
        return Mode::Geometry;
    }

    static bool isNumericArray(const Array& arr)
    {
        ArrayType t = arr.getType();
        return t == ArrayType::DOUBLE || t == ArrayType::SINGLE ||
               t == ArrayType::INT8   || t == ArrayType::UINT8  ||
               t == ArrayType::INT16  || t == ArrayType::UINT16 ||
               t == ArrayType::INT32  || t == ArrayType::UINT32 ||
               t == ArrayType::INT64  || t == ArrayType::UINT64;
    }

    static bool isMatrix2D(const Array& arr)
    {
        auto d = arr.getDimensions();
        return d.size() == 2;
    }

    static std::size_t nRows(const Array& arr)
    {
        return arr.getDimensions()[0];
    }

    static std::size_t nCols(const Array& arr)
    {
        return arr.getDimensions()[1];
    }

    static bool isVectorShape(const Array& arr)
    {
        if (!isMatrix2D(arr)) return false;
        auto d = arr.getDimensions();
        return d[0] == 1 || d[1] == 1;
    }

    static bool isScalarShape(const Array& arr)
    {
        return numel(arr) == 1;
    }

    static bool isLocalPointMatrix(const Array& arr)
    {
        if (!isMatrix2D(arr)) return false;
        if (arr.getType() != ArrayType::DOUBLE) return false;
        std::size_t r = nRows(arr), c = nCols(arr);
        return r >= 3 && (c == 2 || c == 3);
    }

    static bool isNormalVectorForDim(const Array& arr, std::size_t dim)
    {
        if (!isMatrix2D(arr)) return false;
        if (arr.getType() != ArrayType::DOUBLE) return false;
        if (dim == 3) return numel(arr) == 3 && isVectorShape(arr);
        if (dim == 2) return false;
        return false;
    }

    static bool isPossibleNVert(const Array& arr)
    {
        if (!isNumericArray(arr)) return false;
        if (!isVectorShape(arr) && !isScalarShape(arr)) return false;
        return numel(arr) >= 1;
    }

    static bool isBatchNormals(const Array& arr, std::size_t nPoly)
    {
        if (!isMatrix2D(arr)) return false;
        if (arr.getType() != ArrayType::DOUBLE) return false;
        return nRows(arr) == nPoly && nCols(arr) == 3;
    }

    void localGeometry(ArgumentList& outputs, ArgumentList& inputs)
    {
        validateDoubleMatrix(this, factory, inputs[0], "mxPolygonGeometry:input", "points");
        auto dims = inputs[0].getDimensions();
        std::size_t nPts = dims[0], dim = dims[1];
        requireArgs(this, factory, (dim == 2 || dim == 3) && nPts >= 3,
                    "mxPolygonGeometry:input", "points must be N x 2 or N x 3 with N >= 3.");

        std::vector<double> P = typedArrayToVector(inputs[0]);

        std::vector<double> userNormalBuf;
        const double* userNormal = nullptr;
        if (inputs.size() >= 2 && !isCharArray(inputs[1])) {
            requireArgs(this, factory,
                        inputs[1].getType() == ArrayType::DOUBLE &&
                        isVectorShape(inputs[1]) && numel(inputs[1]) == 3,
                        "mxPolygonGeometry:input", "normal must be a real 3-vector.");
            userNormalBuf = typedArrayToVector(inputs[1]);
            userNormal = userNormalBuf.data();
        }

        Mode mode = Mode::Geometry;
        if (inputs.size() >= 2 && isCharArray(inputs[1])) mode = parseMode(inputs[1]);
        if (inputs.size() >= 3)                           mode = parseMode(inputs[2]);

        double area = 0.0, c[3] = {0,0,0}, normal[3] = {0,0,0};
        try {
            if (dim == 2) polygeom::areaCentroidNormal2D(P.data(), nPts, area, c);
            else          polygeom::areaCentroidNormal3D(P.data(), nPts, userNormal, area, c, normal);
        } catch (const std::exception& e) { throwErr("mxPolygonGeometry:error", e.what()); }

        switch (mode) {
            case Mode::Area:
                requireArgs(this, factory, outputs.size() <= 1,
                            "mxPolygonGeometry:output", "Area mode returns one output.");
                outputs[0] = factory.createScalar<double>(area);
                return;

            case Mode::Centroid: {
                requireArgs(this, factory, outputs.size() <= 1,
                            "mxPolygonGeometry:output", "Centroid mode returns one output.");
                std::vector<double> cvec(c, c + dim);
                outputs[0] = vectorToArray(factory, cvec, 1, dim);
                return;
            }

            case Mode::Normal:
                requireArgs(this, factory, outputs.size() <= 1,
                            "mxPolygonGeometry:output", "Normal mode returns one output.");
                requireArgs(this, factory, dim == 3,
                            "mxPolygonGeometry:input", "Normal mode is only defined for N x 3 polygons.");
                outputs[0] = vectorToArray(factory, std::vector<double>(normal, normal + 3), 1, 3);
                return;

            case Mode::Geometry:
                requireArgs(this, factory, outputs.size() <= 3,
                            "mxPolygonGeometry:output",
                            "Geometry mode returns up to three outputs: [area, centroid, normal].");
                outputs[0] = factory.createScalar<double>(area);
                if (outputs.size() >= 2) {
                    std::vector<double> cvec(c, c + dim);
                    outputs[1] = vectorToArray(factory, cvec, 1, dim);
                }
                if (outputs.size() >= 3) {
                    if (dim == 3)
                        outputs[2] = vectorToArray(factory, std::vector<double>(normal, normal + 3), 1, 3);
                    else
                        outputs[2] = factory.createArray<double>({0, 0});
                }
                return;
        }
    }

    void batchGeometry(ArgumentList& outputs, ArgumentList& inputs)
    {
        std::vector<double> Pbuf, nVertBuf;
        polygeom::BatchInput in;
        try {
            in = parseBatchInput(this, factory, inputs[0], inputs[1], Pbuf, nVertBuf);
        } catch (const std::exception& e) { throwErr("mxPolygonGeometry:input", e.what()); }

        std::vector<double> normalsBuf;
        const double* normalsOrNull = nullptr;
        if (inputs.size() >= 3 && !isCharArray(inputs[2])) {
            requireArgs(this, factory, in.dim == 3,
                        "mxPolygonGeometry:input",
                        "Batch normals can only be supplied for 3D polygons.");
            requireArgs(this, factory,
                        inputs[2].getType() == ArrayType::DOUBLE &&
                        getRows(inputs[2]) == in.nPoly && getCols(inputs[2]) == 3,
                        "mxPolygonGeometry:input", "normals must have size nPoly x 3.");
            normalsBuf = typedArrayToVector(inputs[2]);
            normalsOrNull = normalsBuf.data();
        }

        Mode mode = Mode::Geometry;
        if (inputs.size() >= 3 && isCharArray(inputs[2])) mode = parseMode(inputs[2]);
        if (inputs.size() >= 4)                           mode = parseMode(inputs[3]);

        switch (mode) {
            case Mode::Area: {
                requireArgs(this, factory, outputs.size() <= 1,
                            "mxPolygonGeometry:output", "Area mode returns one output.");
                std::vector<double> area;
                try { polygeom::polygonAreaBatch(in, normalsOrNull, area); }
                catch (const std::exception& e) { throwErr("mxPolygonGeometry:error", e.what()); }
                outputs[0] = vectorToArray(factory, area, in.nPoly, 1);
                return;
            }

            case Mode::Centroid: {
                requireArgs(this, factory, outputs.size() <= 1,
                            "mxPolygonGeometry:output", "Centroid mode returns one output.");
                std::vector<double> centroid;
                try { polygeom::polygonCentroidBatch(in, normalsOrNull, centroid); }
                catch (const std::exception& e) { throwErr("mxPolygonGeometry:error", e.what()); }
                outputs[0] = vectorToArray(factory, centroid, in.nPoly, in.dim);
                return;
            }

            case Mode::Normal: {
                requireArgs(this, factory, outputs.size() <= 1,
                            "mxPolygonGeometry:output", "Normal mode returns one output.");
                requireArgs(this, factory, in.dim == 3,
                            "mxPolygonGeometry:input", "Normal mode is only defined for 3D polygons.");
                std::vector<double> normal;
                try { polygeom::polygonNormalBatch(in, normalsOrNull, normal); }
                catch (const std::exception& e) { throwErr("mxPolygonGeometry:error", e.what()); }
                outputs[0] = vectorToArray(factory, normal, in.nPoly, 3);
                return;
            }

            case Mode::Geometry: {
                requireArgs(this, factory, outputs.size() <= 3,
                            "mxPolygonGeometry:output",
                            "Geometry mode returns up to three outputs: [area, centroid, normal].");
                std::vector<double> area, centroid, normal;
                try { polygeom::polygonGeometryBatch(in, normalsOrNull, area, centroid, normal); }
                catch (const std::exception& e) { throwErr("mxPolygonGeometry:error", e.what()); }

                outputs[0] = vectorToArray(factory, area, in.nPoly, 1);
                if (outputs.size() >= 2)
                    outputs[1] = vectorToArray(factory, centroid, in.nPoly, in.dim);
                if (outputs.size() >= 3) {
                    if (in.dim == 3)
                        outputs[2] = vectorToArray(factory, normal, in.nPoly, 3);
                    else
                        outputs[2] = factory.createArray<double>({0, 0});
                }
                return;
            }
        }
    }

     enum class DispatchKind { Local, Batch };

    DispatchKind dispatch(ArgumentList& inputs)
    {
        requireArgs(this, factory, inputs.size() >= 1 && inputs.size() <= 4,
                    "mxPolygonGeometry:input",
                    "Usage: [A,C,N] = mxPolygonGeometry(P [,normal] [,mode]) or "
                    "[A,C,N] = mxPolygonGeometry(Pflat,nVert [,normals] [,mode]).");

        requireArgs(this, factory, isLocalPointMatrix(inputs[0]),
                    "mxPolygonGeometry:input",
                    "First input must be a real double matrix of size N x 2 or N x 3 with N >= 3.");

        const std::size_t dim = nCols(inputs[0]);

        if (inputs.size() == 1) {
            return DispatchKind::Local;
        }

        if (inputs.size() == 2) {
            if (isCharArray(inputs[1])) {
                return DispatchKind::Local;
            }

            if (isNormalVectorForDim(inputs[1], dim)) {
                return DispatchKind::Local;
            }

            if (isPossibleNVert(inputs[1])) {
                return DispatchKind::Batch;
            }

            throwErr("mxPolygonGeometry:input",
                     "Second argument must be either a mode string, a 3-vector normal, or an nVert vector.");
        }

        if (inputs.size() == 3) {
            if (isCharArray(inputs[1])) {
                throwErr("mxPolygonGeometry:input",
                         "Invalid syntax. If the second argument is a mode string, no third argument is allowed.");
            }

            /*
             * Important: decide batch syntax before local-normal syntax.
             * nVert can be scalar, or can itself contain three entries, e.g. [3 4 3].
             * In those cases it can look like a local 3-vector normal unless the
             * third argument is used to disambiguate the call.
             */
            if (isPossibleNVert(inputs[1])) {
                if (isCharArray(inputs[2])) {
                    return DispatchKind::Batch;
                }

                requireArgs(this, factory, dim == 3,
                            "mxPolygonGeometry:input",
                            "Batch normals are only valid for 3D polygons.");

                requireArgs(this, factory, isBatchNormals(inputs[2], numel(inputs[1])),
                            "mxPolygonGeometry:input",
                            "Batch normals must have size nPoly x 3, where nPoly = numel(nVert).");

                return DispatchKind::Batch;
            }

            if (isNormalVectorForDim(inputs[1], dim)) {
                requireArgs(this, factory, isCharArray(inputs[2]),
                            "mxPolygonGeometry:input",
                            "For local syntax with a normal, the third argument must be a mode string.");
                return DispatchKind::Local;
            }

            throwErr("mxPolygonGeometry:input",
                     "Second argument must be either nVert or a 3-vector normal.");
        }

        requireArgs(this, factory, inputs.size() == 4,
                    "mxPolygonGeometry:input", "Too many input arguments.");

        requireArgs(this, factory, !isCharArray(inputs[1]),
                    "mxPolygonGeometry:input",
                    "For 4 inputs, the second argument must be nVert.");

        requireArgs(this, factory, isPossibleNVert(inputs[1]),
                    "mxPolygonGeometry:input",
                    "For 4 inputs, the second argument must be nVert.");

        requireArgs(this, factory, !isCharArray(inputs[2]),
                    "mxPolygonGeometry:input",
                    "For 4 inputs, the third argument must be the normals array.");

        requireArgs(this, factory, dim == 3,
                    "mxPolygonGeometry:input",
                    "Batch normals are only valid for 3D polygons.");

        requireArgs(this, factory, isBatchNormals(inputs[2], numel(inputs[1])),
                    "mxPolygonGeometry:input",
                    "Batch normals must have size nPoly x 3, where nPoly = numel(nVert).");

        requireArgs(this, factory, isCharArray(inputs[3]),
                    "mxPolygonGeometry:input",
                    "For 4 inputs, the fourth argument must be a mode string.");

        return DispatchKind::Batch;
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        DispatchKind kind = dispatch(inputs);
        if (kind == DispatchKind::Local) localGeometry(outputs, inputs);
        else                             batchGeometry(outputs, inputs);
    }
};