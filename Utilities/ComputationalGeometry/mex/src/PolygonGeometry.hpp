#ifndef POLYGON_GEOMETRY_HPP
#define POLYGON_GEOMETRY_HPP

#include "mex.h"
#include <vector>
#include <cstddef>

namespace polygeom {

struct BatchInput {
    const double* P;
    const mwSize* Pdims;
    mwSize nPts;
    mwSize dim;
    const double* nVert;
    mwSize nPoly;
};

void require(bool cond, const char* id, const char* msg);
bool isIntegerValued(double x, double tol = 1e-12);

void rotationFromNormal(const double* n, double* R);

void orderCCW2D(const double* pts, mwSize n, std::vector<mwSize>& perm);
void orderCCW3D(const double* pts, mwSize n, const double* userNormalOrNull, std::vector<mwSize>& perm, double* unitNormalOut);

void areaCentroidNormal2D(const double* pts, mwSize n, double& area, double* centroid);
void areaCentroidNormal3D(const double* pts, mwSize n, double& area, double* centroid, double* unitNormal);

double polygonAreaLocal(const mxArray* points);
mxArray* polygonCentroidLocal(const mxArray* points);
mxArray* polygonNormalLocal(const mxArray* points);
mxArray* orderPointsLocal(int nrhs, const mxArray* prhs[]);

BatchInput parseBatchInput(const mxArray* Pflat, const mxArray* nVert);
void polygonAreaBatch(const BatchInput& in, std::vector<double>& area);
void polygonCentroidBatch(const BatchInput& in, std::vector<double>& centroid);
void polygonNormalBatch(const BatchInput& in, std::vector<double>& normal);
void polygonGeometryBatch(const BatchInput& in,
                          std::vector<double>& area,
                          std::vector<double>& centroid,
                          std::vector<double>& normal);
void orderPointsBatch(const BatchInput& in,
                      std::vector<double>& Pccw,
                      std::vector<double>& perm);

}

#endif
