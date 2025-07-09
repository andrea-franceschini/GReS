#include "mex.h"

void computeQuad4Basis(double* N, const double* coords, mwSize np) {
    for (mwSize i = 0; i < np; ++i) {
        double xi  = coords[i];           // coords(i,0)
        double eta = coords[i + np];      // coords(i,1)

        // Bilinear shape functions for Quad4
        N[i + np * 0] = 0.25 * (1 - xi) * (1 - eta);  // Node 1
        N[i + np * 1] = 0.25 * (1 + xi) * (1 - eta);  // Node 2
        N[i + np * 2] = 0.25 * (1 + xi) * (1 + eta);  // Node 3
        N[i + np * 3] = 0.25 * (1 - xi) * (1 + eta);  // Node 4
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1)
        mexErrMsgTxt("One input required: n×2 coordinate array");

    if (nlhs != 1)
        mexErrMsgTxt("One output required: shape matrix N (n×4)");

    const mxArray *coordList = prhs[0];
    if (!mxIsDouble(coordList) || mxIsComplex(coordList) || mxGetN(coordList) != 2)
        mexErrMsgTxt("Input must be a real double n×2 array");

    mwSize np = mxGetM(coordList);           // number of points
    const double* coords = mxGetPr(coordList); // coords stored column-major

    // Create output array: np × 4
    mwSize dims[2] = {np, 4};
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double* N = mxGetPr(plhs[0]);

    computeQuad4Basis(N, coords, np);
}
