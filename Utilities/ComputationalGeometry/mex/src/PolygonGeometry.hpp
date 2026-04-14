#pragma once

#include "mex.h"
#include <vector>

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

void rotationFromNormal(const double* nIn, double* R);

void orderCCW2D(const double* pts, mwSize n, std::vector<mwSize>& perm);
void orderCCW3D(const double* pts, mwSize n,
                const double* userNormalOrNull,
                std::vector<mwSize>& perm,
                double* unitNormalOut = nullptr);

void areaCentroidNormal2D(const double* pts, mwSize n, double& area, double* centroid);
void areaCentroidNormal3D(const double* pts, mwSize n,
                          const double* userNormalOrNull,
                          double& area, double* centroid, double* unitNormal);

double polygonAreaLocal(const mxArray* points, const mxArray* normalOrNull = nullptr);
mxArray* polygonCentroidLocal(const mxArray* points, const mxArray* normalOrNull = nullptr);
mxArray* polygonNormalLocal(const mxArray* points, const mxArray* normalOrNull = nullptr);
mxArray* orderPointsLocal(int nrhs, const mxArray* prhs[]);

BatchInput parseBatchInput(const mxArray* Pflat, const mxArray* nVert);
const double* parseBatchNormals(const mxArray* normals, mwSize nPoly);

void polygonAreaBatch(const BatchInput& in, const double* normalsOrNull, std::vector<double>& area);
void polygonCentroidBatch(const BatchInput& in, const double* normalsOrNull, std::vector<double>& centroid);
void polygonNormalBatch(const BatchInput& in, const double* normalsOrNull, std::vector<double>& normal);
void polygonGeometryBatch(const BatchInput& in, const double* normalsOrNull,
                          std::vector<double>& area,
                          std::vector<double>& centroid,
                          std::vector<double>& normal);

void orderPointsBatch(const BatchInput& in,
                      const double* normalsOrNull,
                      std::vector<double>& Pccw,
                      std::vector<double>& permOut);

} // namespace polygeom
