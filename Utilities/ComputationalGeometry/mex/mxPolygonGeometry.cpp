#include "mex.h"
#include "src/PolygonGeometry.hpp"

#include <algorithm>
#include <string>
#include <cstring>
#include <vector>

using namespace polygeom;

namespace {

enum class Mode {
    Geometry,
    Area,
    Centroid,
    Normal
};

Mode parseMode(const mxArray* arr) {
    require(mxIsChar(arr), "mxPolygonGeometry:input", "Optional mode must be a character vector or string scalar.");
    char* cstr = mxArrayToString(arr);
    require(cstr != nullptr, "mxPolygonGeometry:input", "Failed to parse mode string.");
    std::string s(cstr);
    mxFree(cstr);

    if (s == "geometry") return Mode::Geometry;
    if (s == "area")     return Mode::Area;
    if (s == "centroid") return Mode::Centroid;
    if (s == "normal")   return Mode::Normal;

    require(false, "mxPolygonGeometry:input", "Mode must be one of: 'geometry', 'area', 'centroid', 'normal'.");
    return Mode::Geometry;
}

bool isLikelyBatchCall(int nrhs, const mxArray* const prhs[]) {
    if (nrhs < 2) return false;
    if (mxIsChar(prhs[1])) return false;
    bool secondIsVector = (mxGetM(prhs[1]) == 1 || mxGetN(prhs[1]) == 1);
    bool firstLooksMatrix = (mxGetN(prhs[0]) == 2 || mxGetN(prhs[0]) == 3);
    bool secondHasMany = (mxGetNumberOfElements(prhs[1]) > 1);
    return firstLooksMatrix && secondIsVector && secondHasMany;
}

} // anonymous namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nrhs >= 1 && nrhs <= 4, "mxPolygonGeometry:input",
            "Usage: [A,C,N] = mxPolygonGeometry(P [,normal] [,mode]) or [A,C,N] = mxPolygonGeometry(Pflat,nVert [,normals] [,mode]).");

    Mode mode = Mode::Geometry;
    bool batchCall = isLikelyBatchCall(nrhs, prhs);

    if (!batchCall) {
        const mxArray* points = prhs[0];
        const mxArray* normalOrNull = nullptr;
        if (nrhs >= 2 && !mxIsChar(prhs[1])) normalOrNull = prhs[1];
        if (nrhs >= 2 && mxIsChar(prhs[1])) mode = parseMode(prhs[1]);
        if (nrhs >= 3) mode = parseMode(prhs[2]);

        switch (mode) {
            case Mode::Area:
                require(nlhs <= 1, "mxPolygonGeometry:output", "Area mode returns one output.");
                plhs[0] = mxCreateDoubleScalar(polygonAreaLocal(points, normalOrNull));
                return;

            case Mode::Centroid:
                require(nlhs <= 1, "mxPolygonGeometry:output", "Centroid mode returns one output.");
                plhs[0] = polygonCentroidLocal(points, normalOrNull);
                return;

            case Mode::Normal:
                require(nlhs <= 1, "mxPolygonGeometry:output", "Normal mode returns one output.");
                plhs[0] = polygonNormalLocal(points, normalOrNull);
                return;

            case Mode::Geometry: {
                require(nlhs <= 3, "mxPolygonGeometry:output", "Geometry mode returns up to three outputs: [area, centroid, normal].");
                double area = polygonAreaLocal(points, normalOrNull);
                mxArray* centroid = polygonCentroidLocal(points, normalOrNull);
                mwSize dim = mxGetN(points);
                mxArray* normal = nullptr;
                if (dim == 3) normal = polygonNormalLocal(points, normalOrNull);

                plhs[0] = mxCreateDoubleScalar(area);
                if (nlhs >= 2) plhs[1] = centroid; else mxDestroyArray(centroid);
                if (nlhs >= 3) {
                    if (dim == 3) plhs[2] = normal;
                    else plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
                } else if (normal) {
                    mxDestroyArray(normal);
                }
                return;
            }
        }
    }

    BatchInput in = parseBatchInput(prhs[0], prhs[1]);
    const double* normalsOrNull = nullptr;

    if (nrhs >= 3 && !mxIsChar(prhs[2])) {
        require(in.dim == 3, "mxPolygonGeometry:input",
                "Batch normals can only be supplied for 3D polygons.");
        normalsOrNull = parseBatchNormals(prhs[2], in.nPoly);
    }

    if (nrhs >= 3 && mxIsChar(prhs[2])) mode = parseMode(prhs[2]);
    if (nrhs >= 4) mode = parseMode(prhs[3]);

    switch (mode) {
        case Mode::Area: {
            require(nlhs <= 1, "mxPolygonGeometry:output", "Area mode returns one output.");
            std::vector<double> area;
            polygonAreaBatch(in, normalsOrNull, area);
            plhs[0] = mxCreateDoubleMatrix(in.nPoly, 1, mxREAL);
            std::copy(area.begin(), area.end(), mxGetPr(plhs[0]));
            return;
        }

        case Mode::Centroid: {
            require(nlhs <= 1, "mxPolygonGeometry:output", "Centroid mode returns one output.");
            std::vector<double> centroid;
            polygonCentroidBatch(in, normalsOrNull, centroid);
            plhs[0] = mxCreateDoubleMatrix(in.nPoly, in.dim, mxREAL);
            std::copy(centroid.begin(), centroid.end(), mxGetPr(plhs[0]));
            return;
        }

        case Mode::Normal: {
            require(nlhs <= 1, "mxPolygonGeometry:output", "Normal mode returns one output.");
            require(in.dim == 3, "mxPolygonGeometry:input", "Normal mode is only defined for 3D polygons.");
            std::vector<double> normal;
            polygonNormalBatch(in, normalsOrNull, normal);
            plhs[0] = mxCreateDoubleMatrix(in.nPoly, 3, mxREAL);
            std::copy(normal.begin(), normal.end(), mxGetPr(plhs[0]));
            return;
        }

        case Mode::Geometry: {
            require(nlhs <= 3, "mxPolygonGeometry:output", "Geometry mode returns up to three outputs: [area, centroid, normal].");
            std::vector<double> area, centroid, normal;
            polygonGeometryBatch(in, normalsOrNull, area, centroid, normal);

            plhs[0] = mxCreateDoubleMatrix(in.nPoly, 1, mxREAL);
            std::copy(area.begin(), area.end(), mxGetPr(plhs[0]));

            if (nlhs >= 2) {
                plhs[1] = mxCreateDoubleMatrix(in.nPoly, in.dim, mxREAL);
                std::copy(centroid.begin(), centroid.end(), mxGetPr(plhs[1]));
            }
            if (nlhs >= 3) {
                if (in.dim == 3) {
                    plhs[2] = mxCreateDoubleMatrix(in.nPoly, 3, mxREAL);
                    std::copy(normal.begin(), normal.end(), mxGetPr(plhs[2]));
                } else {
                    plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
                }
            }
            return;
        }
    }
}
