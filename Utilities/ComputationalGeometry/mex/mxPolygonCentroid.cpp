#include "mex.h"
#include "src/PolygonGeometry.hpp"
#include <vector>
using namespace polygeom;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nlhs <= 1, "mxPolygonCentroid:output", "One output only.");
    if (nrhs == 1) {
        plhs[0] = polygonCentroidLocal(prhs[0]);
        return;
    }
    require(nrhs == 2, "mxPolygonCentroid:input", "Usage: C = mxPolygonCentroid(points) or C = mxPolygonCentroid(Pflat,nVert).");
    BatchInput in = parseBatchInput(prhs[0], prhs[1]);
    std::vector<double> centroid;
    polygonCentroidBatch(in, nullptr, centroid);
    plhs[0] = mxCreateDoubleMatrix(in.nPoly, in.dim, mxREAL);
    std::copy(centroid.begin(), centroid.end(), mxGetPr(plhs[0]));
}
