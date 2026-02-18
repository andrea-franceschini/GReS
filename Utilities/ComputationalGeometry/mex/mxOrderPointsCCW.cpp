#include "mex.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>


struct AngleIndex {
  double angle;
  mwIndex index;
};

inline double dot3(const double* a, const double* b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void cross3(const double* a, const double* b, double* c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

inline double norm3(const double* a) {
  return std::sqrt(dot3(a, a));
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
  /* --- Input checks --- */
  if (nrhs != 2)
    mexErrMsgIdAndTxt("orderPointsCCW:nrhs",
                      "Two inputs required: points, normal.");

  if (nlhs != 1)
    mexErrMsgIdAndTxt("orderPointsCCW:nlhs",
                      "One output required.");

  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    mexErrMsgIdAndTxt("orderPointsCCW:points",
                      "points must be a real double array.");

  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
    mexErrMsgIdAndTxt("orderPointsCCW:normal",
                      "normal must be a real double array.");

  const mwSize* dims = mxGetDimensions(prhs[0]);
  mwSize N = dims[0];
  mwSize dim = dims[1];

  if (dim != 3)
    mexErrMsgIdAndTxt("orderPointsCCW:points",
                      "points must be N x 3.");

  if (mxGetNumberOfElements(prhs[1]) != 3)
    mexErrMsgIdAndTxt("orderPointsCCW:normal",
                      "normal must have 3 elements.");

  const double* points = mxGetPr(prhs[0]);
  const double* normal = mxGetPr(prhs[1]);

  /* --- Compute centroid --- */
  double centroid[3] = {0.0, 0.0, 0.0};

  for (mwSize i = 0; i < N; ++i) {
    centroid[0] += points[i];
    centroid[1] += points[i + N];
    centroid[2] += points[i + 2*N];
  }
  centroid[0] /= N;
  centroid[1] /= N;
  centroid[2] /= N;

  /* --- Reference vector v0 (from point 1) --- */
  double v0[3];
  v0[0] = centroid[0] - points[0];
  v0[1] = centroid[1] - points[N];
  v0[2] = centroid[2] - points[2*N];

  double nrm = norm3(v0);
  if (nrm == 0.0)
    mexErrMsgIdAndTxt("orderPointsCCW:degenerate",
                      "First point coincides with centroid.");

  v0[0] /= nrm;
  v0[1] /= nrm;
  v0[2] /= nrm;

  /* --- Compute angles --- */
  std::vector<AngleIndex> angles(N);
  angles[0] = {0.0, 0};

  for (mwSize i = 1; i < N; ++i) {
    double v[3];
    v[0] = centroid[0] - points[i];
    v[1] = centroid[1] - points[i + N];
    v[2] = centroid[2] - points[i + 2*N];

    double dotv = dot3(v, v0);

    double crossp[3];
    cross3(v, v0, crossp);

    double detv = dot3(normal, crossp);

    angles[i] = { std::atan2(detv, dotv), i };
  }

  std::sort(angles.begin(), angles.end(),
            [](const AngleIndex& a, const AngleIndex& b) {
              return a.angle < b.angle;
            });

  plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
  double* out = mxGetPr(plhs[0]);

  for (mwSize i = 0; i < N; ++i)
    out[i] = static_cast<double>(angles[i].index + 1);
}
