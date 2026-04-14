#include "mex.h"
#include "src/PolygonGeometry.hpp"
#include <vector>
using namespace polygeom;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nlhs <= 1, "mxPolygonNormal:output", "One output only.");
    if (nrhs == 1) {
        plhs[0] = polygonNormalLocal(prhs[0]);
        return;
    }
    require(nrhs == 2, "mxPolygonNormal:input", "Usage: N = mxPolygonNormal(points3D) or N = mxPolygonNormal(Pflat3D,nVert).");
    BatchInput in = parseBatchInput(prhs[0], prhs[1]);
    std::vector<double> normal;
    polygonNormalBatch(in, nullptr, normal);
    plhs[0] = mxCreateDoubleMatrix(in.nPoly, 3, mxREAL);
    std::copy(normal.begin(), normal.end(), mxGetPr(plhs[0]));
}
