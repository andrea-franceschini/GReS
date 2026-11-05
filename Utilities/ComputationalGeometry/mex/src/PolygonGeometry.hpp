#pragma once

#include <cstddef>
#include <vector>
#include <stdexcept>

namespace polygeom {

using Index = std::size_t;

struct BatchInput {
    const double* P;
    Index nPts;
    Index dim;
    const double* nVert;
    Index nPoly;
};

inline void require(bool cond, const char* id, const char* msg) {
    if (!cond) throw std::runtime_error(std::string(id) + ": " + msg);
}

bool isIntegerValued(double x, double tol = 1e-12);

void rotationFromNormal(const double* nIn, double* R);

void orderCCW2D(const double* pts, Index n, std::vector<Index>& perm);
void orderCCW3D(const double* pts, Index n,
                const double* userNormalOrNull,
                std::vector<Index>& perm,
                double* unitNormalOut = nullptr);

void areaCentroidNormal2D(const double* pts, Index n, double& area, double* centroid);
void areaCentroidNormal3D(const double* pts, Index n,
                          const double* userNormalOrNull,
                          double& area, double* centroid, double* unitNormal);

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