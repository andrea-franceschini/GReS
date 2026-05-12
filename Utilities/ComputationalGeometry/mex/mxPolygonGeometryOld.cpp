#include "mex.hpp"
#include "mexAdapter.hpp"
#include "src/MexHelper.hpp"

#include <vector>
#include <algorithm>
#include <string>

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
                        inputs[1].getType() == ArrayType::DOUBLE && numel(inputs[1]) == 3,
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

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        requireArgs(this, factory, inputs.size() >= 1 && inputs.size() <= 4,
                    "mxPolygonGeometry:input",
                    "Usage: [A,C,N] = mxPolygonGeometry(P [,normal] [,mode]) or "
                    "[A,C,N] = mxPolygonGeometry(Pflat,nVert [,normals] [,mode]).");

        bool batchCall = false;
        if (inputs.size() >= 2) {
            bool secondIsVector = (getRows(inputs[1]) == 1 || getCols(inputs[1]) == 1);
            bool firstLooksMatrix = (getCols(inputs[0]) == 2 || getCols(inputs[0]) == 3);
            bool secondHasMany = (numel(inputs[1]) > 1);
            batchCall = !isCharArray(inputs[1]) && firstLooksMatrix && secondIsVector && secondHasMany;
        }

        if (!batchCall) localGeometry(outputs, inputs);
        else            batchGeometry(outputs, inputs);
    }
};