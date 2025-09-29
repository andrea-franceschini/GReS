#include "mex.h"
#include <cmath>

void computeBasis(double* N, const double* coords, mwSize np) {
    for (mwSize i = 0; i < np; ++i) {
        double xi  = coords[i];           // xi = c(i,0)
        double eta = coords[i + np];      // eta = c(i,1)

        // 1D basis functions
        double b1_xi  = 0.5 * xi * (xi - 1);
        double b2_xi  = 1.0 - xi * xi;
        double b3_xi  = 0.5 * xi * (xi + 1);

        double b1_eta = 0.5 * eta * (eta - 1);
        double b2_eta = 1.0 - eta * eta;
        double b3_eta = 0.5 * eta * (eta + 1);

        // Compute tensor product basis functions
        N[i + np * 0] = b1_xi * b1_eta;  // Node 1
        N[i + np * 1] = b3_xi * b1_eta;  // Node 2
        N[i + np * 2] = b3_xi * b3_eta;  // Node 3
        N[i + np * 3] = b1_xi * b3_eta;  // Node 4
        N[i + np * 4] = b2_xi * b1_eta;  // Node 5
        N[i + np * 5] = b3_xi * b2_eta;  // Node 6
        N[i + np * 6] = b2_xi * b3_eta;  // Node 7
        N[i + np * 7] = b1_xi * b2_eta;  // Node 8
        N[i + np * 8] = b2_xi * b2_eta;  // Node 9
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1)
        mexErrMsgTxt("One input required: coordList (n×2 array)");
    if (nlhs != 1)
        mexErrMsgTxt("One output required: N (n×9 array)");

    const mxArray *coordList = prhs[0];
    if (!mxIsDouble(coordList) || mxIsComplex(coordList) || mxGetN(coordList) != 2)
        mexErrMsgTxt("Input must be a real double n×2 array");

    mwSize np = mxGetM(coordList); // number of points
    const double* coords = mxGetPr(coordList); // column-major

    // Create output array: np × 9
    mwSize dims[2] = {np, 9};
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double* N = mxGetPr(plhs[0]);

    computeBasis(N, coords, np);
}
