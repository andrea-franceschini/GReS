#include "mex.h"
#include "src/PolygonGeometry.hpp"
#include <vector>
using namespace polygeom;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nrhs == 2, "mxPolygonGeometry:input", "Usage: [A,C,N] = mxPolygonGeometry(Pflat,nVert). N is returned only for 3D.");
    BatchInput in = parseBatchInput(prhs[0], prhs[1]);
    std::vector<double> area, centroid, normal;
    polygonGeometryBatch(in, area, centroid, normal);

    require(nlhs >= 2 && nlhs <= 3, "mxPolygonGeometry:output", "Use [A,C] for 2D and [A,C,N] for 3D (or [A,C] in 3D if normal is not needed).");

    plhs[0] = mxCreateDoubleMatrix(in.nPoly, 1, mxREAL);
    std::copy(area.begin(), area.end(), mxGetPr(plhs[0]));

    plhs[1] = mxCreateDoubleMatrix(in.nPoly, in.dim, mxREAL);
    std::copy(centroid.begin(), centroid.end(), mxGetPr(plhs[1]));

    if (nlhs >= 3) {
        require(in.dim == 3, "mxPolygonGeometry:output", "Normal output is only available for 3D input polygons.");
        plhs[2] = mxCreateDoubleMatrix(in.nPoly, 3, mxREAL);
        std::copy(normal.begin(), normal.end(), mxGetPr(plhs[2]));
    }
}
