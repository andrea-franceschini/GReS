#include "PolygonGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <cstddef>

namespace polygeom {

namespace {

inline double dot2(const double* a, const double* b) {
    return a[0]*b[0] + a[1]*b[1];
}
inline double dot3(const double* a, const double* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline void cross3(const double* a, const double* b, double* c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}
inline double norm2(const double* a) {
    return std::sqrt(dot2(a,a));
}
inline double norm3(const double* a) {
    return std::sqrt(dot3(a,a));
}
inline void normalize3(double* a) {
    double n = norm3(a);
    require(n > 1e-15, "PolygonGeometry:degenerate", "Zero or near-zero vector cannot be normalized.");
    a[0] /= n; a[1] /= n; a[2] /= n;
}
inline void getPoint(const double* P, std::size_t nRows, std::size_t dim, std::size_t i, double* out) {
    out[0] = P[i];
    out[1] = P[i + nRows];
    if (dim == 3) out[2] = P[i + 2*nRows];
}

struct AngleIdx {
    double ang;
    std::size_t idx;
};

void computeNewellNormalOrdered(const double* pts, std::size_t n, double* normal) {
    normal[0] = 0.0; normal[1] = 0.0; normal[2] = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t j = (i + 1) % n;
        double pi[3], pj[3];
        getPoint(pts, n, 3, i, pi);
        getPoint(pts, n, 3, j, pj);
        normal[0] += (pi[1] - pj[1]) * (pi[2] + pj[2]);
        normal[1] += (pi[2] - pj[2]) * (pi[0] + pj[0]);
        normal[2] += (pi[0] - pj[0]) * (pi[1] + pj[1]);
    }
    normalize3(normal);
}

bool normalFromAnyTriple(const double* pts, std::size_t n, double* normal) {
    double p0[3], p1[3], p2[3], v1[3], v2[3], cp[3];

    for (std::size_t i = 0; i < n; ++i) {
        getPoint(pts, n, 3, i, p0);
        for (std::size_t j = i + 1; j < n; ++j) {
            getPoint(pts, n, 3, j, p1);
            v1[0] = p1[0] - p0[0];
            v1[1] = p1[1] - p0[1];
            v1[2] = p1[2] - p0[2];
            for (std::size_t k = j + 1; k < n; ++k) {
                getPoint(pts, n, 3, k, p2);
                v2[0] = p2[0] - p0[0];
                v2[1] = p2[1] - p0[1];
                v2[2] = p2[2] - p0[2];
                cross3(v1, v2, cp);
                double nrm = norm3(cp);
                if (nrm > 1e-15) {
                    normal[0] = cp[0] / nrm;
                    normal[1] = cp[1] / nrm;
                    normal[2] = cp[2] / nrm;
                    return true;
                }
            }
        }
    }
    return false;
}

void orthonormalBasisFromNormal(const double* n, double* e1, double* e2) {
    double tmp[3];
    if (std::fabs(n[0]) < 0.9) { tmp[0]=1.0; tmp[1]=0.0; tmp[2]=0.0; }
    else                       { tmp[0]=0.0; tmp[1]=1.0; tmp[2]=0.0; }

    double proj = dot3(tmp, n);
    e1[0] = tmp[0] - proj*n[0];
    e1[1] = tmp[1] - proj*n[1];
    e1[2] = tmp[2] - proj*n[2];
    normalize3(e1);
    cross3(n, e1, e2);
    normalize3(e2);
}

void project3Dto2D(const double* pts3, std::size_t n, const double* normal, std::vector<double>& pts2) {
    double e1[3], e2[3];
    orthonormalBasisFromNormal(normal, e1, e2);
    pts2.assign(n * 2, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double p[3]; getPoint(pts3, n, 3, i, p);
        pts2[i] = dot3(p, e1);
        pts2[i + n] = dot3(p, e2);
    }
}

void orderFromProjected2D(const std::vector<double>& pts2, std::size_t n, std::vector<std::size_t>& perm) {
    double c[2] = {0.0, 0.0};
    for (std::size_t i = 0; i < n; ++i) {
        c[0] += pts2[i];
        c[1] += pts2[i + n];
    }
    c[0] /= static_cast<double>(n);
    c[1] /= static_cast<double>(n);

    std::vector<AngleIdx> ang(n);
    for (std::size_t i = 0; i < n; ++i) {
        double x = pts2[i] - c[0];
        double y = pts2[i + n] - c[1];
        ang[i].ang = std::atan2(y, x);
        ang[i].idx = i;
    }
    std::sort(ang.begin(), ang.end(), [](const AngleIdx& a, const AngleIdx& b) { return a.ang < b.ang; });
    perm.resize(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = ang[i].idx;

    double twiceArea = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t ia = perm[i];
        std::size_t ib = perm[(i+1)%n];
        double xa = pts2[ia], ya = pts2[ia + n];
        double xb = pts2[ib], yb = pts2[ib + n];
        twiceArea += xa*yb - xb*ya;
    }
    if (twiceArea < 0.0) std::reverse(perm.begin(), perm.end());
}

void reorderPoints(const double* ptsIn, std::size_t n, std::size_t dim, const std::vector<std::size_t>& perm, std::vector<double>& ptsOut) {
    ptsOut.assign(n * dim, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t src = perm[i];
        ptsOut[i] = ptsIn[src];
        ptsOut[i + n] = ptsIn[src + n];
        if (dim == 3) ptsOut[i + 2*n] = ptsIn[src + 2*n];
    }
}

void areaCentroid2DOrdered(const double* pts, std::size_t n, double& area, double* centroid) {
    if (n == 3) {
        double x1 = pts[0],     y1 = pts[0 + n];
        double x2 = pts[1],     y2 = pts[1 + n];
        double x3 = pts[2],     y3 = pts[2 + n];
        double cross = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
        area = 0.5 * std::fabs(cross);
        require(area > 1e-15, "PolygonGeometry:degenerate", "Zero-area triangle.");
        centroid[0] = (x1 + x2 + x3) / 3.0;
        centroid[1] = (y1 + y2 + y3) / 3.0;
        return;
    }
    if (n == 4) {
        double tri1[6] = {pts[0], pts[1], pts[2], pts[n+0], pts[n+1], pts[n+2]};
        double tri2[6] = {pts[0], pts[2], pts[3], pts[n+0], pts[n+2], pts[n+3]};
        double A1, A2, c1[2], c2[2];
        areaCentroid2DOrdered(tri1, 3, A1, c1);
        areaCentroid2DOrdered(tri2, 3, A2, c2);
        area = A1 + A2;
        centroid[0] = (A1*c1[0] + A2*c2[0]) / area;
        centroid[1] = (A1*c1[1] + A2*c2[1]) / area;
        return;
    }

    double A2 = 0.0, Cx = 0.0, Cy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        std::size_t j = (i + 1) % n;
        double xi = pts[i], yi = pts[i + n];
        double xj = pts[j], yj = pts[j + n];
        double cr = xi*yj - xj*yi;
        A2 += cr;
        Cx += (xi + xj) * cr;
        Cy += (yi + yj) * cr;
    }
    require(std::fabs(A2) > 1e-15, "PolygonGeometry:degenerate", "Zero-area polygon.");
    area = 0.5 * std::fabs(A2);
    centroid[0] = Cx / (3.0 * A2);
    centroid[1] = Cy / (3.0 * A2);
}

void areaCentroid3DOrdered(const double* pts, std::size_t n, double& area, double* centroid, double* unitNormal) {
    if (n == 3) {
        double p0[3], p1[3], p2[3], v1[3], v2[3], cp[3];
        getPoint(pts, n, 3, 0, p0); getPoint(pts, n, 3, 1, p1); getPoint(pts, n, 3, 2, p2);
        v1[0]=p1[0]-p0[0]; v1[1]=p1[1]-p0[1]; v1[2]=p1[2]-p0[2];
        v2[0]=p2[0]-p0[0]; v2[1]=p2[1]-p0[1]; v2[2]=p2[2]-p0[2];
        cross3(v1,v2,cp);
        area = 0.5 * norm3(cp);
        require(area > 1e-15, "PolygonGeometry:degenerate", "Zero-area triangle.");
        unitNormal[0]=cp[0]; unitNormal[1]=cp[1]; unitNormal[2]=cp[2];
        normalize3(unitNormal);
        centroid[0]=(p0[0]+p1[0]+p2[0])/3.0;
        centroid[1]=(p0[1]+p1[1]+p2[1])/3.0;
        centroid[2]=(p0[2]+p1[2]+p2[2])/3.0;
        return;
    }
    if (n == 4) {
        std::vector<double> tri1(9), tri2(9);
        for (int d = 0; d < 3; ++d) {
            tri1[0 + 3*d] = pts[0 + n*d];
            tri1[1 + 3*d] = pts[1 + n*d];
            tri1[2 + 3*d] = pts[2 + n*d];
            tri2[0 + 3*d] = pts[0 + n*d];
            tri2[1 + 3*d] = pts[2 + n*d];
            tri2[2 + 3*d] = pts[3 + n*d];
        }
        double A1,A2,c1[3],c2[3],n1[3],n2[3];
        areaCentroid3DOrdered(tri1.data(),3,A1,c1,n1);
        areaCentroid3DOrdered(tri2.data(),3,A2,c2,n2);
        area = A1 + A2;
        centroid[0] = (A1*c1[0] + A2*c2[0]) / area;
        centroid[1] = (A1*c1[1] + A2*c2[1]) / area;
        centroid[2] = (A1*c1[2] + A2*c2[2]) / area;
        unitNormal[0] = n1[0] + n2[0];
        unitNormal[1] = n1[1] + n2[1];
        unitNormal[2] = n1[2] + n2[2];
        normalize3(unitNormal);
        return;
    }

    computeNewellNormalOrdered(pts, n, unitNormal);
    double p0[3]; getPoint(pts, n, 3, 0, p0);
    double C[3] = {0.0, 0.0, 0.0};
    area = 0.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        double p1[3], p2[3], v1[3], v2[3], cp[3], ct[3];
        getPoint(pts, n, 3, i, p1);
        getPoint(pts, n, 3, i+1, p2);
        v1[0]=p1[0]-p0[0]; v1[1]=p1[1]-p0[1]; v1[2]=p1[2]-p0[2];
        v2[0]=p2[0]-p0[0]; v2[1]=p2[1]-p0[1]; v2[2]=p2[2]-p0[2];
        cross3(v1,v2,cp);
        double a = 0.5 * norm3(cp);
        require(a > 1e-15, "PolygonGeometry:degenerate", "Degenerate triangulation in polygon.");
        ct[0]=(p0[0]+p1[0]+p2[0])/3.0;
        ct[1]=(p0[1]+p1[1]+p2[1])/3.0;
        ct[2]=(p0[2]+p1[2]+p2[2])/3.0;
        C[0] += a*ct[0]; C[1] += a*ct[1]; C[2] += a*ct[2];
        area += a;
    }
    require(area > 1e-15, "PolygonGeometry:degenerate", "Zero-area polygon.");
    centroid[0] = C[0]/area; centroid[1] = C[1]/area; centroid[2] = C[2]/area;
}

inline const double* getBatchNormalPtr(const double* normalsOrNull, std::size_t p, std::size_t nPoly, double* buf) {
    if (!normalsOrNull) return nullptr;
    buf[0] = normalsOrNull[p];
    buf[1] = normalsOrNull[p + nPoly];
    buf[2] = normalsOrNull[p + 2*nPoly];
    return buf;
}

} // anonymous namespace


bool isIntegerValued(double x, double tol) {
    return std::fabs(x - std::round(x)) <= tol;
}

void rotationFromNormal(const double* nIn, double* R) {
    double n[3] = {nIn[0], nIn[1], nIn[2]};
    normalize3(n);
    double e1[3], e2[3];
    orthonormalBasisFromNormal(n, e1, e2);
    for (int i = 0; i < 3; ++i) {
        R[i + 0*3] = n[i];
        R[i + 1*3] = e1[i];
        R[i + 2*3] = e2[i];
    }
}

void orderCCW2D(const double* pts, std::size_t n, std::vector<std::size_t>& perm) {
    std::vector<double> pts2(2*n);
    std::copy(pts, pts + 2*n, pts2.begin());
    orderFromProjected2D(pts2, n, perm);
}

void orderCCW3D(const double* pts, std::size_t n,
                const double* userNormalOrNull,
                std::vector<std::size_t>& perm,
                double* unitNormalOut) {
    double normal[3];
    if (userNormalOrNull) {
        normal[0] = userNormalOrNull[0];
        normal[1] = userNormalOrNull[1];
        normal[2] = userNormalOrNull[2];
        normalize3(normal);
    } else {
        bool ok = normalFromAnyTriple(pts, n, normal);
        require(ok, "PolygonGeometry:degenerate",
                "Could not determine a valid plane normal from the input points.");
    }

    std::vector<double> pts2;
    project3Dto2D(pts, n, normal, pts2);
    orderFromProjected2D(pts2, n, perm);

    if (unitNormalOut) {
        unitNormalOut[0] = normal[0];
        unitNormalOut[1] = normal[1];
        unitNormalOut[2] = normal[2];
    }
}

void areaCentroidNormal2D(const double* pts, std::size_t n, double& area, double* centroid) {
    std::vector<std::size_t> perm;
    orderCCW2D(pts, n, perm);
    std::vector<double> ordered;
    reorderPoints(pts, n, 2, perm, ordered);
    areaCentroid2DOrdered(ordered.data(), n, area, centroid);
}

void areaCentroidNormal3D(const double* pts, std::size_t n,
                          const double* userNormalOrNull,
                          double& area, double* centroid, double* unitNormal) {
    std::vector<std::size_t> perm;
    double normal0[3];
    orderCCW3D(pts, n, userNormalOrNull, perm, normal0);

    std::vector<double> ordered;
    reorderPoints(pts, n, 3, perm, ordered);

    areaCentroid3DOrdered(ordered.data(), n, area, centroid, unitNormal);
    if (dot3(unitNormal, normal0) < 0.0) {
        unitNormal[0] = -unitNormal[0];
        unitNormal[1] = -unitNormal[1];
        unitNormal[2] = -unitNormal[2];
    }
}

void polygonAreaBatch(const BatchInput& in, const double* normalsOrNull, std::vector<double>& area) {
    area.assign(in.nPoly, 0.0);
    if (normalsOrNull) {
        require(in.dim == 3, "PolygonGeometry:input", "A normals array can only be provided for 3D polygons.");
    }

    std::size_t off = 0;
    for (std::size_t p = 0; p < in.nPoly; ++p) {
        std::size_t n = static_cast<std::size_t>(std::llround(in.nVert[p]));
        std::vector<double> poly(n * in.dim);
        for (std::size_t i = 0; i < n; ++i) {
            poly[i] = in.P[off + i];
            poly[i + n] = in.P[off + i + in.nPts];
            if (in.dim == 3) poly[i + 2*n] = in.P[off + i + 2*in.nPts];
        }
        double c[3], normal[3], normalBuf[3];
        if (in.dim == 2) areaCentroidNormal2D(poly.data(), n, area[p], c);
        else {
            const double* userNormal = getBatchNormalPtr(normalsOrNull, p, in.nPoly, normalBuf);
            areaCentroidNormal3D(poly.data(), n, userNormal, area[p], c, normal);
        }
        off += n;
    }
}

void polygonCentroidBatch(const BatchInput& in, const double* normalsOrNull, std::vector<double>& centroid) {
    centroid.assign(in.nPoly * in.dim, 0.0);
    if (normalsOrNull) {
        require(in.dim == 3, "PolygonGeometry:input", "A normals array can only be provided for 3D polygons.");
    }

    std::size_t off = 0;
    for (std::size_t p = 0; p < in.nPoly; ++p) {
        std::size_t n = static_cast<std::size_t>(std::llround(in.nVert[p]));
        std::vector<double> poly(n * in.dim);
        for (std::size_t i = 0; i < n; ++i) {
            poly[i] = in.P[off + i];
            poly[i + n] = in.P[off + i + in.nPts];
            if (in.dim == 3) poly[i + 2*n] = in.P[off + i + 2*in.nPts];
        }
        double area, c[3], normal[3], normalBuf[3];
        if (in.dim == 2) areaCentroidNormal2D(poly.data(), n, area, c);
        else {
            const double* userNormal = getBatchNormalPtr(normalsOrNull, p, in.nPoly, normalBuf);
            areaCentroidNormal3D(poly.data(), n, userNormal, area, c, normal);
        }
        centroid[p] = c[0];
        centroid[p + in.nPoly] = c[1];
        if (in.dim == 3) centroid[p + 2*in.nPoly] = c[2];
        off += n;
    }
}

void polygonNormalBatch(const BatchInput& in, const double* normalsOrNull, std::vector<double>& normal) {
    require(in.dim == 3, "PolygonGeometry:input", "Polygon normals are only defined here for 3D polygons.");
    normal.assign(in.nPoly * 3, 0.0);

    std::size_t off = 0;
    for (std::size_t p = 0; p < in.nPoly; ++p) {
        std::size_t n = static_cast<std::size_t>(std::llround(in.nVert[p]));
        std::vector<double> poly(n * 3);
        for (std::size_t i = 0; i < n; ++i) {
            poly[i] = in.P[off + i];
            poly[i + n] = in.P[off + i + in.nPts];
            poly[i + 2*n] = in.P[off + i + 2*in.nPts];
        }
        double area, c[3], nn[3], normalBuf[3];
        const double* userNormal = getBatchNormalPtr(normalsOrNull, p, in.nPoly, normalBuf);
        areaCentroidNormal3D(poly.data(), n, userNormal, area, c, nn);
        normal[p] = nn[0];
        normal[p + in.nPoly] = nn[1];
        normal[p + 2*in.nPoly] = nn[2];
        off += n;
    }
}

void polygonGeometryBatch(const BatchInput& in,
                          const double* normalsOrNull,
                          std::vector<double>& area,
                          std::vector<double>& centroid,
                          std::vector<double>& normal) {
    area.assign(in.nPoly, 0.0);
    centroid.assign(in.nPoly * in.dim, 0.0);
    if (in.dim == 3) normal.assign(in.nPoly * 3, 0.0); else normal.clear();

    if (normalsOrNull) {
        require(in.dim == 3, "PolygonGeometry:input", "A normals array can only be provided for 3D polygons.");
    }

    std::size_t off = 0;
    for (std::size_t p = 0; p < in.nPoly; ++p) {
        std::size_t n = static_cast<std::size_t>(std::llround(in.nVert[p]));
        std::vector<double> poly(n * in.dim);
        for (std::size_t i = 0; i < n; ++i) {
            poly[i] = in.P[off + i];
            poly[i + n] = in.P[off + i + in.nPts];
            if (in.dim == 3) poly[i + 2*n] = in.P[off + i + 2*in.nPts];
        }
        double c[3] = {0.0,0.0,0.0};
        double nn[3] = {0.0,0.0,0.0};
        double normalBuf[3];
        if (in.dim == 2) areaCentroidNormal2D(poly.data(), n, area[p], c);
        else {
            const double* userNormal = getBatchNormalPtr(normalsOrNull, p, in.nPoly, normalBuf);
            areaCentroidNormal3D(poly.data(), n, userNormal, area[p], c, nn);
        }
        centroid[p] = c[0];
        centroid[p + in.nPoly] = c[1];
        if (in.dim == 3) centroid[p + 2*in.nPoly] = c[2];
        if (in.dim == 3) {
            normal[p] = nn[0];
            normal[p + in.nPoly] = nn[1];
            normal[p + 2*in.nPoly] = nn[2];
        }
        off += n;
    }
}

void orderPointsBatch(const BatchInput& in,
                      const double* normalsOrNull,
                      std::vector<double>& Pccw,
                      std::vector<double>& permOut) {
    Pccw.assign(in.nPts * in.dim, 0.0);
    permOut.assign(in.nPts, 0.0);

    if (normalsOrNull) {
        require(in.dim == 3, "PolygonGeometry:input",
                "A normals array can only be provided for 3D polygons.");
    }

    std::size_t off = 0;
    for (std::size_t p = 0; p < in.nPoly; ++p) {
        std::size_t n = static_cast<std::size_t>(std::llround(in.nVert[p]));
        std::vector<double> poly(n * in.dim);
        for (std::size_t i = 0; i < n; ++i) {
            poly[i] = in.P[off + i];
            poly[i + n] = in.P[off + i + in.nPts];
            if (in.dim == 3) poly[i + 2*n] = in.P[off + i + 2*in.nPts];
        }
        std::vector<std::size_t> perm;
        if (in.dim == 2) {
            orderCCW2D(poly.data(), n, perm);
        } else {
            double normalBuf[3];
            const double* userNormal = getBatchNormalPtr(normalsOrNull, p, in.nPoly, normalBuf);
            orderCCW3D(poly.data(), n, userNormal, perm, nullptr);
        }
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t src = perm[i];
            Pccw[off + i] = poly[src];
            Pccw[off + i + in.nPts] = poly[src + n];
            if (in.dim == 3) Pccw[off + i + 2*in.nPts] = poly[src + 2*n];
            permOut[off + i] = static_cast<double>(src + 1);
        }
        off += n;
    }
}

} // namespace polygeom