#include "mex.h"
#include <array>
#include <cmath>

// cross product
inline std::array<double,3> cross(const std::array<double,3>& a,
                                  const std::array<double,3>& b) {
    return { a[1]*b[2]-a[2]*b[1],
             a[2]*b[0]-a[0]*b[2],
             a[0]*b[1]-a[1]*b[0] };
}

// dot product
inline double dot(const std::array<double,3>& a, const std::array<double,3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// norm
inline double norm(const std::array<double,3>& a) {
    return std::sqrt(dot(a,a));
}

// normalize
inline std::array<double,3> normalize(const std::array<double,3>& a) {
    double n = norm(a);
    if(n<1e-16) mexErrMsgTxt("Zero-length vector cannot be normalized.");
    return { a[0]/n, a[1]/n, a[2]/n };
}

// determinant of 3x3 matrix given as columns
inline double det3x3(const std::array<std::array<double,3>,3>& R) {
    const auto& a = R[0]; const auto& b = R[1]; const auto& c = R[2];
    return a[0]*(b[1]*c[2]-b[2]*c[1])
         - a[1]*(b[0]*c[2]-b[2]*c[0])
         + a[2]*(b[0]*c[1]-b[1]*c[0]);
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if(nrhs!=1) mexErrMsgTxt("Usage: R = computeRot_mex(n)");
    if(nlhs!=1) mexErrMsgTxt("One output required.");

    if(mxGetNumberOfElements(prhs[0])!=3 || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("Input n must be a 3-element real vector.");

    double* n_ptr = mxGetPr(prhs[0]);
    std::array<double,3> n = { n_ptr[0], n_ptr[1], n_ptr[2] };
    n = normalize(n);

    // pick a vector not parallel to n
    std::array<double,3> tmp;
    if(std::fabs(n[0])<0.9) tmp = {1.0,0.0,0.0};
    else tmp = {0.0,1.0,0.0};

    // first tangent: orthogonalize tmp against n
    double ndot = dot(tmp,n);
    std::array<double,3> m1 = { tmp[0]-ndot*n[0],
                                tmp[1]-ndot*n[1],
                                tmp[2]-ndot*n[2] };
    m1 = normalize(m1);

    // second tangent: orthogonal to both
    std::array<double,3> m2 = cross(n,m1);
    m2 = normalize(m2);

    // assemble rotation matrix: columns = n, m1, m2
    plhs[0] = mxCreateDoubleMatrix(3,3,mxREAL);
    double* R_ptr = mxGetPr(plhs[0]);

    std::array<std::array<double,3>,3> Rcols = { n, m1, m2 };

    // enforce right-handed system (det=+1)
    double detR = det3x3(Rcols);
    if(detR<0.0){
        m1[0] = -m1[0]; m1[1] = -m1[1]; m1[2] = -m1[2];
        Rcols = { n, m1, m2 };
        detR = det3x3(Rcols);
    }

    if(std::fabs(detR-1.0)>1e-12) mexErrMsgTxt("Rotation matrix is not orthogonal to machine precision.");

    // fill MATLAB column-major array
    for(int j=0;j<3;++j)
        for(int i=0;i<3;++i)
            R_ptr[i + j*3] = Rcols[j][i];
}
