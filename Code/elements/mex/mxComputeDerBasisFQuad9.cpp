#include "mex.h"
#include <cmath>

void computeDerivatives(double* dN, const double* coords, mwSize np) {
    for (mwSize i = 0; i < np; ++i) {
        double xi  = coords[i];           // c(i,0)
        double eta = coords[i + np];      // c(i,1)

        // Basis functions
        double b1_xi  = 0.5 * xi * (xi - 1);
        double b2_xi  = 1.0 - xi * xi;
        double b3_xi  = 0.5 * xi * (xi + 1);
        double gb1_xi = 0.5 * (2 * xi - 1);
        double gb2_xi = -2.0 * xi;
        double gb3_xi = 0.5 * (2 * xi + 1);

        double b1_eta  = 0.5 * eta * (eta - 1);
        double b2_eta  = 1.0 - eta * eta;
        double b3_eta  = 0.5 * eta * (eta + 1);
        double gb1_eta = 0.5 * (2 * eta - 1);
        double gb2_eta = -2.0 * eta;
        double gb3_eta = 0.5 * (2 * eta + 1);

        // Indexing: dN(row, col, i) => dN[row + 2 * col + 18 * i]
        // where row = 0 (d/dxi), 1 (d/deta); col = 0..8

        // Node 1
        dN[0 + 2 * 0 + 18 * i] = gb1_xi * b1_eta;
        dN[1 + 2 * 0 + 18 * i] = b1_xi  * gb1_eta;

        // Node 2
        dN[0 + 2 * 1 + 18 * i] = gb3_xi * b1_eta;
        dN[1 + 2 * 1 + 18 * i] = b3_xi  * gb1_eta;

        // Node 3
        dN[0 + 2 * 2 + 18 * i] = gb3_xi * b3_eta;
        dN[1 + 2 * 2 + 18 * i] = b3_xi  * gb3_eta;

        // Node 4
        dN[0 + 2 * 3 + 18 * i] = gb1_xi * b3_eta;
        dN[1 + 2 * 3 + 18 * i] = b1_xi  * gb3_eta;

        // Node 5
        dN[0 + 2 * 4 + 18 * i] = gb2_xi * b1_eta;
        dN[1 + 2 * 4 + 18 * i] = b2_xi  * gb1_eta;

        // Node 6
        dN[0 + 2 * 5 + 18 * i] = gb3_xi * b2_eta;
        dN[1 + 2 * 5 + 18 * i] = b3_xi  * gb2_eta;

        // Node 7
        dN[0 + 2 * 6 + 18 * i] = gb2_xi * b3_eta;
        dN[1 + 2 * 6 + 18 * i] = b2_xi  * gb3_eta;

        // Node 8
        dN[0 + 2 * 7 + 18 * i] = gb1_xi * b2_eta;
        dN[1 + 2 * 7 + 18 * i] = b1_xi  * gb2_eta;

        // Node 9
        dN[0 + 2 * 8 + 18 * i] = gb2_xi * b2_eta;
        dN[1 + 2 * 8 + 18 * i] = b2_xi  * gb2_eta;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1)
        mexErrMsgTxt("One input required: list (n×2 array)");
    if (nlhs != 1)
        mexErrMsgTxt("One output required: dN (2×9×n array)");

    const mxArray *list = prhs[0];
    if (!mxIsDouble(list) || mxIsComplex(list) || mxGetN(list) != 2)
        mexErrMsgTxt("Input must be a real double n×2 array");

    mwSize np = mxGetM(list);
    const double* coords = mxGetPr(list);  // coords stored column-major

    // Output array: 2×9×np
    mwSize dims[3] = {2, 9, np};
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double* dN = mxGetPr(plhs[0]);

    computeDerivatives(dN, coords, np);
}
