#include "mex.h"
#include "src/PolygonGeometry.hpp"
#include <vector>
using namespace polygeom;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nlhs <= 1, "mxPolygonArea:output", "One output only.");
    if (nrhs == 1) {
        plhs[0] = mxCreateDoubleScalar(polygonAreaLocal(prhs[0]));
        return;
    }
    require(nrhs == 2, "mxPolygonArea:input", "Usage: A = mxPolygonArea(points) or A = mxPolygonArea(Pflat,nVert).");
    BatchInput in = parseBatchInput(prhs[0], prhs[1]);
    std::vector<double> area;
    polygonAreaBatch(in, nullptr, area);
    plhs[0] = mxCreateDoubleMatrix(in.nPoly, 1, mxREAL);
    std::copy(area.begin(), area.end(), mxGetPr(plhs[0]));
}
