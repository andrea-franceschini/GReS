#include "mex.h"
#include <cmath>
#include <vector>

// Compute determinant of 3x3 matrix (column-major)
double det3x3(const double* M) {
  return M[0]*(M[4]*M[8] - M[5]*M[7])
       - M[1]*(M[3]*M[8] - M[5]*M[6])
       + M[2]*(M[3]*M[7] - M[4]*M[6]);
}

// Compute inverse of 3x3 matrix (column-major)
void inv3x3(const double* M, double* invM) {
  double det = det3x3(M);
  if (fabs(det) < 1e-15) {
    mexErrMsgIdAndTxt("mexFunction:singularMatrix", "Jacobian matrix is singular.");
  }
  double invDet = 1.0 / det;

  invM[0] =  (M[4]*M[8] - M[5]*M[7]) * invDet;
  invM[1] = -(M[1]*M[8] - M[2]*M[7]) * invDet;
  invM[2] =  (M[1]*M[5] - M[2]*M[4]) * invDet;
  invM[3] = -(M[3]*M[8] - M[5]*M[6]) * invDet;
  invM[4] =  (M[0]*M[8] - M[2]*M[6]) * invDet;
  invM[5] = -(M[0]*M[5] - M[2]*M[3]) * invDet;
  invM[6] =  (M[3]*M[7] - M[4]*M[6]) * invDet;
  invM[7] = -(M[0]*M[7] - M[1]*M[6]) * invDet;
  invM[8] =  (M[0]*M[4] - M[1]*M[3]) * invDet;
}

// Multiply 3xN by Nx3 = 3x3
void mul3xN_Nx3(const double* A, const double* B, double* out, int N) {
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      double sum = 0.0;
      for (int k = 0; k < N; ++k) {
        sum += A[row + 3*k] * B[k + N*col];
      }
      out[row + 3*col] = sum;
    }
  }
}

// Multiply 3x3 by 3xN = 3xN
void mul3x3_3xN(const double* A, const double* B, double* out, int N) {
  for (int col = 0; col < N; ++col) {
    for (int row = 0; row < 3; ++row) {
      double sum = 0.0;
      for (int k = 0; k < 3; ++k) {
        sum += A[row + 3*k] * B[k + 3*col];
      }
      out[row + 3*col] = sum;
    }
  }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
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
    mexErrMsgIdAndTxt("mexFunction:invalidInput",
                      "J1 must be of size 3 x Nnodes x nGauss");
  }

  int Nnodes = dimsJ1[1];
  int nGauss = (ndimsJ1 >= 3) ? dimsJ1[2] : 1;

  const mwSize* dimsCoords = mxGetDimensions(coords_mx);
  if (mxGetNumberOfDimensions(coords_mx) != 2 ||
      dimsCoords[0] != Nnodes || dimsCoords[1] != 3) {
    mexErrMsgIdAndTxt("mexFunction:invalidInput",
                      "coords must be Nnodes x 3");
  }

  const mwSize* dimsW = mxGetDimensions(weights_mx);
  if (!((dimsW[0] == nGauss && dimsW[1] == 1) ||
        (dimsW[1] == nGauss && dimsW[0] == 1))) {
    mexErrMsgIdAndTxt("mexFunction:invalidInput",
                      "weights must be a vector of length nGauss");
  }

  const double* J1 = mxGetPr(J1_mx);
  const double* coords = mxGetPr(coords_mx);
  const double* weights = mxGetPr(weights_mx);

  // Create output arrays
  mwSize dimsN[3] = {3, (mwSize)Nnodes, (mwSize)nGauss};
  plhs[0] = mxCreateNumericArray(3, dimsN, mxDOUBLE_CLASS, mxREAL);
  double* N = mxGetPr(plhs[0]);

  mwSize dimsDet[2] = {static_cast<mwSize>(nGauss), 1};

  plhs[1] = mxCreateNumericArray(2, dimsDet, mxDOUBLE_CLASS, mxREAL);
  double* detJ = mxGetPr(plhs[1]);

  // Temporary buffers
  double J[9];
  double invJ[9];
  std::vector<double> Ntemp(3 * Nnodes);

  for (int gp = 0; gp < nGauss; ++gp) {
    const double* J1gp = J1 + gp * 3 * Nnodes;

    // Compute Jacobian J = J1(:,:,gp) * coords
    mul3xN_Nx3(J1gp, coords, J, Nnodes);

    // Determinant × weight
    detJ[gp] = det3x3(J) * weights[gp];

    // Inverse
    inv3x3(J, invJ);

    // Compute N(:,:,gp) = invJ * J1(:,:,gp)
    mul3x3_3xN(invJ, J1gp, Ntemp.data(), Nnodes);

    // Store result
    for (int i = 0; i < 3 * Nnodes; ++i) {
      N[gp * 3 * Nnodes + i] = Ntemp[i];
    }
  }
}
