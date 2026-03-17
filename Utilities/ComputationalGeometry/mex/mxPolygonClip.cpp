#include "mex.h"
#include "polygonClip.hpp"

// Convert MATLAB Nx2 array -> Polygon
Polygon mxToPolygon(const mxArray* mx){
    mwSize N=mxGetM(mx);
    const double* P=mxGetPr(mx);
    Polygon poly(N);
    for(mwSize i=0;i<N;++i) poly[i]={P[i],P[i+N]};
    return poly;
}

// Convert Polygon -> MATLAB Nx2 array
mxArray* polygonToMx(const Polygon& poly){
    mwSize N=poly.size();
    mxArray* mx=mxCreateDoubleMatrix(N,2,mxREAL);
    double* P=mxGetPr(mx);
    for(mwSize i=0;i<N;++i){ P[i]=poly[i][0]; P[i+N]=poly[i][1]; }
    return mx;
}

void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[]){
    if(nrhs!=2) mexErrMsgTxt("Usage: [polyOut,isValid] = polygonClip_mex(poly, clipPoly)");
    if(nlhs!=2) mexErrMsgTxt("Two outputs required.");

    for(int i=0;i<2;++i){
        if(!mxIsDouble(prhs[i])||mxIsComplex(prhs[i])||mxGetN(prhs[i])!=2)
            mexErrMsgTxt("Inputs must be Nx2 real double matrices.");
    }

    Polygon poly    = mxToPolygon(prhs[0]);
    Polygon clipper = mxToPolygon(prhs[1]);

    bool valid=false;
    if(poly.size()>=3 && clipper.size()>=3)
        valid=clipPolygon(poly,clipper);

    plhs[0]=polygonToMx(valid?poly:Polygon{});   // clipped polygon
    plhs[1]=mxCreateLogicalScalar(valid);       // validity flag
}
