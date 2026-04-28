#include "mex.h"
#include <math.h>

#if defined(_MSC_VER)
#define RESTRICT __restrict
#define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#define FORCE_INLINE static inline __attribute__((always_inline))
#else
#define RESTRICT
#define FORCE_INLINE static inline
#endif

// Unified determinant and inverse computation. Removes redundant cofactor calculations.
FORCE_INLINE void compute_inv_det(const double *RESTRICT M, double *RESTRICT invM, double *RESTRICT det_out) {
    invM[0] = M[4]*M[8] - M[5]*M[7];
    invM[1] = M[2]*M[7] - M[1]*M[8];
    invM[2] = M[1]*M[5] - M[2]*M[4];
    invM[3] = M[5]*M[6] - M[3]*M[8];
    invM[4] = M[0]*M[8] - M[2]*M[6];
    invM[5] = M[2]*M[3] - M[0]*M[5];
    invM[6] = M[3]*M[7] - M[4]*M[6];
    invM[7] = M[1]*M[6] - M[0]*M[7];
    invM[8] = M[0]*M[4] - M[1]*M[3];

    double det = M[0]*invM[0] + M[3]*invM[1] + M[6]*invM[2];
    
    if (fabs(det) < 1e-15) {
        mexErrMsgIdAndTxt("mexFunction:singularMatrix", "Jacobian matrix is singular.");
    }

    double invDet = 1.0 / det;
    for(int i = 0; i < 9; ++i) {
        invM[i] *= invDet;
    }
    *det_out = det;
}

// Core loop. If N is passed as a literal constant, inner loops are fully unrolled at compile time.
FORCE_INLINE void process_gauss_points(int N, int nGauss, 
                                       const double *RESTRICT J1, 
                                       const double *RESTRICT coords, 
                                       const double *RESTRICT weights, 
                                       double *RESTRICT N_out, 
                                       double *RESTRICT detJ) 
{
    double J[9];
    double invJ[9];

    for (int gp = 0; gp < nGauss; ++gp) {
        const double *RESTRICT J1gp = J1 + gp * 3 * N;
        double *RESTRICT N_target = N_out + gp * 3 * N;

        // Multiply 3xN by Nx3 -> 3x3
        for (int col = 0; col < 3; ++col) {
            for (int row = 0; row < 3; ++row) {
                double sum = 0.0;
                for (int k = 0; k < N; ++k) {
                    sum += J1gp[row + 3*k] * coords[k + N*col];
                }
                J[row + 3*col] = sum;
            }
        }

        // Determinant and Inverse
        double det;
        compute_inv_det(J, invJ, &det);
        detJ[gp] = det * weights[gp];

        // Multiply 3x3 by 3xN -> 3xN (Zero-copy directly to output buffer)
        for (int col = 0; col < N; ++col) {
            for (int row = 0; row < 3; ++row) {
                double sum = 0.0;
                for (int k = 0; k < 3; ++k) {
                    sum += invJ[row + 3*k] * J1gp[k + 3*col];
                }
                N_target[row + 3*col] = sum;
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("mexFunction:invalidNumInputs",
                          "3 inputs required: J1(3xNnodesxN), coords(Nnodesx3), weights(N)");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("mexFunction:invalidNumOutputs",
                          "2 outputs required: N (3xNnodesxN), detJ (Nx1)");
    }

    const mxArray* J1_mx = prhs[0];
    const mxArray* coords_mx = prhs[1];
    const mxArray* weights_mx = prhs[2];

    const mwSize* dimsJ1 = mxGetDimensions(J1_mx);
    int ndimsJ1 = mxGetNumberOfDimensions(J1_mx);
    if (ndimsJ1 < 2 || dimsJ1[0] != 3) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "J1 must be of size 3 x Nnodes x nGauss");
    }

    int Nnodes = dimsJ1[1];
    int nGauss = (ndimsJ1 >= 3) ? dimsJ1[2] : 1;

    const mwSize* dimsCoords = mxGetDimensions(coords_mx);
    if (mxGetNumberOfDimensions(coords_mx) != 2 || dimsCoords[0] != Nnodes || dimsCoords[1] != 3) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "coords must be Nnodes x 3");
    }

    const mwSize* dimsW = mxGetDimensions(weights_mx);
    if (!((dimsW[0] == nGauss && dimsW[1] == 1) || (dimsW[1] == nGauss && dimsW[0] == 1))) {
        mexErrMsgIdAndTxt("mexFunction:invalidInput", "weights must be a vector of length nGauss");
    }

    const double *RESTRICT J1 = mxGetPr(J1_mx);
    const double *RESTRICT coords = mxGetPr(coords_mx);
    const double *RESTRICT weights = mxGetPr(weights_mx);

    mwSize dimsN[3] = {3, (mwSize)Nnodes, (mwSize)nGauss};
    plhs[0] = mxCreateNumericArray(3, dimsN, mxDOUBLE_CLASS, mxREAL);
    double *RESTRICT N_out = mxGetPr(plhs[0]);

    mwSize dimsDet[2] = {(mwSize)(nGauss), 1};
    plhs[1] = mxCreateNumericArray(2, dimsDet, mxDOUBLE_CLASS, mxREAL);
    double *RESTRICT detJ = mxGetPr(plhs[1]);

    // Explicit compile-time constant branching for aggressive unrolling
    if (Nnodes == 4) {         // Tetra
        process_gauss_points(4, nGauss, J1, coords, weights, N_out, detJ);
    } else if (Nnodes == 8) {  // Hexa8
        process_gauss_points(8, nGauss, J1, coords, weights, N_out, detJ);
    } else if (Nnodes == 20) { // Hexa20
        process_gauss_points(20, nGauss, J1, coords, weights, N_out, detJ);
    } else if (Nnodes == 27) { // Hexa27
        process_gauss_points(27, nGauss, J1, coords, weights, N_out, detJ);
    } else {                   // Arbitrary element fallback
        process_gauss_points(Nnodes, nGauss, J1, coords, weights, N_out, detJ);
    }
}