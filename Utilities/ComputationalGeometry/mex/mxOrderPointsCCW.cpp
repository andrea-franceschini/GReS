#include "mex.h"
#include "src/PolygonGeometry.hpp"

#include <algorithm>
#include <vector>

using namespace polygeom;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nlhs <= 2, "mxOrderPointsCCW:output",
            "At most two outputs: [idx] or [Pccw, perm].");
    require(nrhs >= 1 && nrhs <= 3, "mxOrderPointsCCW:input",
            "Usage: idx = mxOrderPointsCCW(P [,normal]) or [Pccw,perm] = mxOrderPointsCCW(Pflat,nVert [,normals]).");

    bool batchCall = false;
    if (nrhs >= 2) {
        bool secondIsVector = (mxGetM(prhs[1]) == 1 || mxGetN(prhs[1]) == 1);
        bool firstLooksMatrix = (mxGetN(prhs[0]) == 2 || mxGetN(prhs[0]) == 3);
        bool secondHasMany = (mxGetNumberOfElements(prhs[1]) > 1);
        batchCall = firstLooksMatrix && secondIsVector && secondHasMany;
    }

    if (!batchCall) {
        require(nlhs <= 1, "mxOrderPointsCCW:output", "Local call returns only idx.");
        plhs[0] = orderPointsLocal(nrhs, prhs);
        return;
    }

    BatchInput in = parseBatchInput(prhs[0], prhs[1]);

    const double* normalsOrNull = nullptr;
    if (nrhs == 3) {
        require(in.dim == 3, "mxOrderPointsCCW:input",
                "Batch normals can only be supplied for 3D polygons.");
        normalsOrNull = parseBatchNormals(prhs[2], in.nPoly);
    }

    std::vector<double> Pccw, perm;
    orderPointsBatch(in, normalsOrNull, Pccw, perm);

    plhs[0] = mxCreateDoubleMatrix(in.nPts, in.dim, mxREAL);
    std::copy(Pccw.begin(), Pccw.end(), mxGetPr(plhs[0]));

    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(in.nPts, 1, mxREAL);
        std::copy(perm.begin(), perm.end(), mxGetPr(plhs[1]));
    }
}
