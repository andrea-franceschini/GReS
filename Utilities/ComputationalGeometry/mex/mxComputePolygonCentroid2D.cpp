#include "mex.h"
#include <cmath>

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
  /* ---- checks ---- */
  if (nrhs != 1)
    mexErrMsgTxt("Usage: centroid = computePolygonCentroid2D_CCW_mex(points)");

  if (nlhs != 1)
    mexErrMsgTxt("One output required.");

  const mxArray* Pmx = prhs[0];

  if (!mxIsDouble(Pmx) || mxIsComplex(Pmx))
    mexErrMsgTxt("points must be real double.");

  mwSize N = mxGetM(Pmx);
  mwSize dim = mxGetN(Pmx);

  if (dim != 2)
    mexErrMsgTxt("points must be N x 2.");

  if (N < 3)
    mexErrMsgTxt("At least 3 points required.");

  const double* P = mxGetPr(Pmx);

  /* ---- centroid computation ---- */
  double Cx = 0.0;
  double Cy = 0.0;
  double A  = 0.0;   // signed area * 2

  for (mwSize i = 0; i < N; ++i)
  {
    mwSize j = (i + 1) % N;

    double xi = P[i];
    double yi = P[i + N];

    double xj = P[j];
    double yj = P[j + N];

    double cross = xi*yj - xj*yi;

    A  += cross;
    Cx += (xi + xj) * cross;
    Cy += (yi + yj) * cross;
  }

  if (A == 0.0)
    mexErrMsgTxt("Zero-area polygon.");

  A *= 0.5;

  Cx /= (6.0 * A);
  Cy /= (6.0 * A);

  /* ---- output ---- */
  plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
  double* out = mxGetPr(plhs[0]);

  out[0] = Cx;
  out[1] = Cy;
}
