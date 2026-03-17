#include "mex.h"
#include <cmath>

inline double dot3(const double* a, const double* b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void cross3(const double* a, const double* b, double* c)
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

inline double norm3(const double* a)
{
  return std::sqrt(dot3(a, a));
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
  // checks
  if (nrhs != 1)
    mexErrMsgTxt("Usage: centroid = computePolygonCentroidCCW_mex(points)");

  if (nlhs != 1)
    mexErrMsgTxt("One output required.");

  const mxArray* Pmx = prhs[0];

  if (!mxIsDouble(Pmx) || mxIsComplex(Pmx))
    mexErrMsgTxt("points must be real double.");

  mwSize N = mxGetM(Pmx);
  if (mxGetN(Pmx) != 3)
    mexErrMsgTxt("points must be N x 3.");

  if (N < 3)
    mexErrMsgTxt("At least 3 points required.");

  const double* P = mxGetPr(Pmx);

  // reference point
  double P0[3] = { P[0], P[N], P[2*N] };

  double C[3] = {0.0, 0.0, 0.0};
  double Atot = 0.0;

  // fan triangulation
  for (mwSize i = 1; i < N-1; ++i)
  {
    double v1[3] = {
      P[i]     - P0[0],
      P[i + N] - P0[1],
      P[i + 2*N] - P0[2]
    };

    double v2[3] = {
      P[i+1]     - P0[0],
      P[i+1 + N] - P0[1],
      P[i+1 + 2*N] - P0[2]
    };

    double cp[3];
    cross3(v1, v2, cp);
    double A = 0.5 * norm3(cp);

    double Ct[3] = {
      (P0[0] + P[i] + P[i+1]) / 3.0,
      (P0[1] + P[i + N] + P[i+1 + N]) / 3.0,
      (P0[2] + P[i + 2*N] + P[i+1 + 2*N]) / 3.0
    };

    C[0] += A * Ct[0];
    C[1] += A * Ct[1];
    C[2] += A * Ct[2];
    Atot += A;
  }

  if (Atot == 0.0)
    mexErrMsgTxt("Zero-area polygon.");

  C[0] /= Atot;
  C[1] /= Atot;
  C[2] /= Atot;

  /* ---- output ---- */
  plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);
  double* out = mxGetPr(plhs[0]);
  out[0] = C[0];
  out[1] = C[1];
  out[2] = C[2];
}
