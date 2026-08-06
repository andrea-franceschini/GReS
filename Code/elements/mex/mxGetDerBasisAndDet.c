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

/* Compute the inverse and determinant of a column-major 3x3 matrix. */
FORCE_INLINE void compute_inv_det(const double *RESTRICT M,
                                  double *RESTRICT invM,
                                  double *RESTRICT det_out)
{
  invM[0] = M[4]*M[8] - M[5]*M[7];
  invM[1] = M[2]*M[7] - M[1]*M[8];
  invM[2] = M[1]*M[5] - M[2]*M[4];
  invM[3] = M[5]*M[6] - M[3]*M[8];
  invM[4] = M[0]*M[8] - M[2]*M[6];
  invM[5] = M[2]*M[3] - M[0]*M[5];
  invM[6] = M[3]*M[7] - M[4]*M[6];
  invM[7] = M[1]*M[6] - M[0]*M[7];
  invM[8] = M[0]*M[4] - M[1]*M[3];

  const double det = M[0]*invM[0] + M[3]*invM[1] + M[6]*invM[2];

  if (fabs(det) < 1e-15) {
    mexErrMsgIdAndTxt("mxGetDerBasisAndDet:singularMatrix",
                      "A Jacobian matrix is singular.");
  }

  const double invDet = 1.0 / det;
  for (int i = 0; i < 9; ++i) {
    invM[i] *= invDet;
  }

  *det_out = det;
}

/*
* J1:      3 x N x nGauss, shared by all cells
* coords:  N x 3 x nCells
* N_out:   3 x N x nGauss x nCells
* detJ:    nGauss x nCells
*/
FORCE_INLINE void process_cells(int N,
                                mwSize nGauss,
                                mwSize nCells,
                                const double *RESTRICT J1,
                                const double *RESTRICT coords,
                                const double *RESTRICT weights,
                                double *RESTRICT N_out,
                                double *RESTRICT detJ)
{
  const mwSize j1GaussStride = (mwSize)3 * (mwSize)N;
  const mwSize coordsCellStride = (mwSize)N * 3;
  const mwSize outputGaussStride = j1GaussStride;
  const mwSize outputCellStride = outputGaussStride * nGauss;

  for (mwSize cell = 0; cell < nCells; ++cell) {
    const double *RESTRICT coordsCell = coords + cell * coordsCellStride;
    double *RESTRICT outputCell = N_out + cell * outputCellStride;
    double *RESTRICT detCell = detJ + cell * nGauss;

    for (mwSize gp = 0; gp < nGauss; ++gp) {
      const double *RESTRICT J1gp = J1 + gp * j1GaussStride;
      double *RESTRICT N_target = outputCell + gp * outputGaussStride;
      double J[9];
      double invJ[9];

      /* J = J1gp * coordsCell: (3 x N)(N x 3) -> 3 x 3. */
      for (int col = 0; col < 3; ++col) {
        for (int row = 0; row < 3; ++row) {
          double sum = 0.0;
          for (int k = 0; k < N; ++k) {
            sum += J1gp[row + 3*k] * coordsCell[k + N*col];
          }
          J[row + 3*col] = sum;
        }
      }

      double det;
      compute_inv_det(J, invJ, &det);
      detCell[gp] = det * weights[gp];

      /* N_target = invJ * J1gp: (3 x 3)(3 x N) -> 3 x N. */
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
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs != 3) {
    mexErrMsgIdAndTxt("mxGetDerBasisAndDet:invalidNumInputs",
                      "Three inputs required: J1(3xNnodesxnGauss), "
                      "coords(Nnodesx3xnCells), weights(nGauss).");
  }
  if (nlhs != 2) {
    mexErrMsgIdAndTxt("mxGetDerBasisAndDet:invalidNumOutputs",
                      "Two outputs required: N(3xNnodesxnGaussxnCells), "
                      "detJ(nGaussxnCells).");
  }

  const mxArray *J1_mx = prhs[0];
  const mxArray *coords_mx = prhs[1];
  const mxArray *weights_mx = prhs[2];

  if (!mxIsDouble(J1_mx) || mxIsComplex(J1_mx) ||
      !mxIsDouble(coords_mx) || mxIsComplex(coords_mx) ||
      !mxIsDouble(weights_mx) || mxIsComplex(weights_mx)) {
    mexErrMsgIdAndTxt("mxGetDerBasisAndDet:invalidType",
                      "All inputs must be real double arrays.");
  }

  const mwSize ndimsJ1 = mxGetNumberOfDimensions(J1_mx);
  const mwSize *dimsJ1 = mxGetDimensions(J1_mx);
  if (ndimsJ1 > 3 || dimsJ1[0] != 3) {
    mexErrMsgIdAndTxt("mxGetDerBasisAndDet:invalidJ1",
                      "J1 must have size 3 x Nnodes x nGauss.");
  }

  const mwSize Nnodes = dimsJ1[1];
  const mwSize nGauss = (ndimsJ1 >= 3) ? dimsJ1[2] : 1;

  const mwSize ndimsCoords = mxGetNumberOfDimensions(coords_mx);
  const mwSize *dimsCoords = mxGetDimensions(coords_mx);
  if (ndimsCoords > 3 || dimsCoords[0] != Nnodes || dimsCoords[1] != 3) {
    mexErrMsgIdAndTxt("mxGetDerBasisAndDet:invalidCoords",
                      "coords must have size Nnodes x 3 x nCells.");
  }
  const mwSize nCells = (ndimsCoords >= 3) ? dimsCoords[2] : 1;

  if (mxGetNumberOfElements(weights_mx) != nGauss) {
    mexErrMsgIdAndTxt("mxGetDerBasisAndDet:invalidWeights",
                      "weights must contain exactly nGauss elements.");
  }

  const double *RESTRICT J1 = mxGetPr(J1_mx);
  const double *RESTRICT coords = mxGetPr(coords_mx);
  const double *RESTRICT weights = mxGetPr(weights_mx);

  const mwSize dimsN[4] = {3, Nnodes, nGauss, nCells};
  plhs[0] = mxCreateNumericArray(4, dimsN, mxDOUBLE_CLASS, mxREAL);
  double *RESTRICT N_out = mxGetPr(plhs[0]);

  const mwSize dimsDet[2] = {nGauss, nCells};
  plhs[1] = mxCreateNumericArray(2, dimsDet, mxDOUBLE_CLASS, mxREAL);
  double *RESTRICT detJ = mxGetPr(plhs[1]);

  /* Preserve constant element-size branches for compiler unrolling. */
  if (Nnodes == 4) {
    process_cells(4, nGauss, nCells, J1, coords, weights, N_out, detJ);
  } else if (Nnodes == 8) {
    process_cells(8, nGauss, nCells, J1, coords, weights, N_out, detJ);
  } else if (Nnodes == 20) {
    process_cells(20, nGauss, nCells, J1, coords, weights, N_out, detJ);
  } else if (Nnodes == 27) {
    process_cells(27, nGauss, nCells, J1, coords, weights, N_out, detJ);
  } else {
    process_cells((int)Nnodes, nGauss, nCells,
                  J1, coords, weights, N_out, detJ);
  }
}