#include "mex.h"

void computeQuad4Derivatives(double* dN, const double* coords, mwSize np) {
    for (mwSize i = 0; i < np; ++i) {
        double xi  = coords[i];           // coords(i,0)
        double eta = coords[i + np];      // coords(i,1)

        // Derivatives dN_dxi and dN_deta for each node
        double dN_dxi[4] = {
            -0.25 * (1 - eta),   // Node 1
             0.25 * (1 - eta),   // Node 2
             0.25 * (1 + eta),   // Node 3
            -0.25 * (1 + eta)    // Node 4
        };

        double dN_deta[4] = {
            -0.25 * (1 - xi),    // Node 1
            -0.25 * (1 + xi),    // Node 2
             0.25 * (1 + xi),    // Node 3
             0.25 * (1 - xi)     // Node 4
        };

        for (mwSize j = 0; j < 4; ++j) {
            dN[0 + 2 * j + 8 * i] = dN_dxi[j];   // d/dxi
            dN[1 + 2 * j + 8 * i] = dN_deta[j];  // d/deta
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1)
        mexErrMsgTxt("One input required: n×2 coordinates");

    if (nlhs != 1)
        mexErrMsgTxt("One output required: 2×4×n array");

    const mxArray *list = prhs[0];
    if (!mxIsDouble(list) || mxIsComplex(list) || mxGetN(list) != 2)
        mexErrMsgTxt("Input must be a real double n×2 array");

    mwSize np = mxGetM(list);
    const double* coords = mxGetPr(list);

    mwSize dims[3] = {2, 4, np};  // 2 (xi, eta) × 4 nodes × np points
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* dN = mxGetPr(plhs[0]);

    computeQuad4Derivatives(dN, coords, np);
}
