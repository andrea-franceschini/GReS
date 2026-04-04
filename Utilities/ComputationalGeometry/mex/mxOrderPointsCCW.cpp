#include "mex.h"
#include "src/PolygonGeometry.hpp"
#include <vector>
using namespace polygeom;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nlhs <= 2, "mxOrderPointsCCW:output", "At most two outputs: [idx] or [Pccw,perm].");

    bool batchCall = (nrhs == 2 && mxGetN(prhs[0]) >= 2 && mxGetM(prhs[1]) * mxGetN(prhs[1]) > 1);
    if (!batchCall) {
        require(nlhs <= 1, "mxOrderPointsCCW:output", "Local call returns only idx.");
        plhs[0] = orderPointsLocal(nrhs, prhs);
        return;
    }

    BatchInput in = parseBatchInput(prhs[0], prhs[1]);
    std::vector<double> Pccw, perm;
    orderPointsBatch(in, Pccw, perm);
    plhs[0] = mxCreateDoubleMatrix(in.nPts, in.dim, mxREAL);
    std::copy(Pccw.begin(), Pccw.end(), mxGetPr(plhs[0]));
    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(in.nPts, 1, mxREAL);
        std::copy(perm.begin(), perm.end(), mxGetPr(plhs[1]));
    }
}
