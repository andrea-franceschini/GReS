#include "mex.h"
#include "src/PolygonGeometry.hpp"
using namespace polygeom;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nrhs == 1, "mxComputeRotationMat:input", "Usage: R = mxComputeRotationMat(normal)." );
    require(nlhs <= 1, "mxComputeRotationMat:output", "One output only.");
    require(mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) && mxGetNumberOfElements(prhs[0]) == 3,
            "mxComputeRotationMat:input", "normal must be a real 3-vector.");
    const double* n = mxGetPr(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(3,3,mxREAL);
    rotationFromNormal(n, mxGetPr(plhs[0]));
}
