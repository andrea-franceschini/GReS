#include "mex.h"

extern "C" {
#include "gmshc.h"
}

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

// Gmsh owns the allocation routines for C API output arrays.
template<class T> struct GmshBuffer {
    T *data = nullptr;
    size_t size = 0;

    ~GmshBuffer() { gmshFree(data); }
};

struct FragmentMap {
    GmshBuffer<int *> entities;
    size_t *sizes = nullptr;

    ~FragmentMap()
    {
        for (size_t i = 0; i < entities.size; ++i)
            gmshFree(entities.data[i]);
        gmshFree(sizes);
    }
};

void checkGmsh(int code, const char *call)
{
    if (code != 0)
        throw std::runtime_error(std::string(call) + " failed (error " +
                                 std::to_string(code) + ").");
}

int addPolygon(const mxArray *polygon)
{
    const size_t n = mxGetM(polygon);
    const double *xyz = mxGetPr(polygon);
    std::vector<int> points(n), lines(n);
    int ierr = 0;
    for (size_t i = 0; i < n; ++i) {
        points[i] = gmshModelOccAddPoint(xyz[i], xyz[n + i], xyz[2 * n + i], 0, -1, &ierr);
        checkGmsh(ierr, "gmshModelOccAddPoint");
    }
    for (size_t i = 0; i < n; ++i) {
        lines[i] = gmshModelOccAddLine(points[i], points[(i + 1) % n], -1, &ierr);
        checkGmsh(ierr, "gmshModelOccAddLine");
    }

    const int loop = gmshModelOccAddCurveLoop(lines.data(), n, -1, &ierr);
    checkGmsh(ierr, "gmshModelOccAddCurveLoop");
    const int surface = gmshModelOccAddPlaneSurface(&loop, 1, -1, &ierr);
    checkGmsh(ierr, "gmshModelOccAddPlaneSurface");
    return surface;
}

struct Mesh {
    std::vector<double> coordinates, tetrahedra;
    std::vector<std::vector<double>> fractures;
};

void appendElements(int type, int entity, const std::unordered_map<size_t, size_t> &nodeRows,
                    std::vector<double> &output)
{
    GmshBuffer<size_t> elements, nodes;
    int ierr = 0;
    gmshModelMeshGetElementsByType(type, &elements.data, &elements.size,
                                   &nodes.data, &nodes.size, entity, 0, 1, &ierr);
    checkGmsh(ierr, "gmshModelMeshGetElementsByType");
    for (size_t i = 0; i < nodes.size; ++i)
        output.push_back(static_cast<double>(nodeRows.at(nodes.data[i])));
}

struct MeshSizing {
    double matrix;
    double fracture;
    double transition;
};

using Point = std::array<double, 3>;

Point subtract(const Point &a, const Point &b)
{
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double dot(const Point &a, const Point &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Point cross(const Point &a, const Point &b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

struct Grading {
    std::vector<std::vector<Point>> polygons;
    std::vector<Point> normals;
    double matrixSize, fractureSize, nearDistance, width;
};

// Project onto the polygon plane. Use the plane distance inside the outline
// and the distance to the closest perimeter segment outside it.
double distanceToPolygon(const Point &point, const std::vector<Point> &polygon,
                         const Point &normal)
{
    const double height = dot(subtract(point, polygon.front()), normal);
    Point projected = point;
    for (size_t axis = 0; axis < 3; ++axis)
        projected[axis] -= height * normal[axis];

    size_t axis = 0;
    for (size_t i = 1; i < 3; ++i)
        if (std::abs(normal[i]) > std::abs(normal[axis]))
            axis = i;
    const size_t u = (axis + 1) % 3, v = (axis + 2) % 3;
    bool inside = false;
    double nearestSquared = INFINITY;
    for (size_t i = 0; i < polygon.size(); ++i) {
        const Point &a = polygon[i], &b = polygon[(i + 1) % polygon.size()];
        const Point edge = subtract(b, a), offset = subtract(point, a);
        const double t = std::max(0., std::min(1., dot(offset, edge) / dot(edge, edge)));
        const Point delta{offset[0] - t * edge[0], offset[1] - t * edge[1],
                          offset[2] - t * edge[2]};
        nearestSquared = std::min(nearestSquared, dot(delta, delta));
        if ((a[v] > projected[v]) != (b[v] > projected[v]) &&
            projected[u] < a[u] + (projected[v] - a[v]) * edge[u] / edge[v])
            inside = !inside;
    }
    return inside ? std::abs(height) : std::sqrt(nearestSquared);
}

double meshSizeAt(int, int, double x, double y, double z, double, void *context)
{
    const auto &grading = *static_cast<const Grading *>(context);
    const Point point{x, y, z};
    double distance = grading.nearDistance + grading.width;
    for (size_t i = 0; i < grading.polygons.size(); ++i)
        distance = std::min(distance, distanceToPolygon(point, grading.polygons[i], grading.normals[i]));

    const double weight = std::max(0., (distance - grading.nearDistance) / grading.width);
    return grading.fractureSize + weight * (grading.matrixSize - grading.fractureSize);
}

void generateMesh(const MeshSizing &sizing, const mxArray *polygons)
{
    const double matrixSize = sizing.matrix > 0 ? sizing.matrix : sizing.fracture;
    const double fractureSize = sizing.fracture > 0 ? sizing.fracture : matrixSize;
    const size_t nFractures = mxGetNumberOfElements(polygons);
    const bool graded = nFractures > 0 && matrixSize != fractureSize;
    Grading grading;
    int ierr = 0;
    if (graded) {
        grading.matrixSize = matrixSize;
        grading.fractureSize = fractureSize;
        grading.nearDistance = 0.5 * std::min(matrixSize, fractureSize);
        grading.width = sizing.transition * std::abs(matrixSize - fractureSize);
        if (!std::isfinite(grading.nearDistance + grading.width) || grading.width <= 0)
            throw std::runtime_error("Invalid mesh grading width.");

        // The input polygons cover the same area as their OCC fragments.
        for (size_t i = 0; i < nFractures; ++i) {
            const mxArray *input = mxGetCell(polygons, i);
            const size_t n = mxGetM(input);
            const double *xyz = mxGetPr(input);
            std::vector<Point> polygon(n);
            for (size_t j = 0; j < n; ++j)
                polygon[j] = {xyz[j], xyz[n + j], xyz[2 * n + j]};

            Point normal{};
            for (size_t j = 1; j + 1 < n; ++j) {
                const Point wedge = cross(subtract(polygon[j], polygon.front()),
                                          subtract(polygon[j + 1], polygon.front()));
                for (size_t axis = 0; axis < 3; ++axis)
                    normal[axis] += wedge[axis];
            }
            const double length = std::sqrt(dot(normal, normal));
            if (!(length > 0) || !std::isfinite(length))
                throw std::runtime_error("Invalid fracture normal for mesh grading.");
            for (double &component : normal)
                component /= length;
            grading.polygons.push_back(std::move(polygon));
            grading.normals.push_back(normal);
        }

        gmshOptionSetNumber("Mesh.MeshSizeExtendFromBoundary", 0, &ierr);
        checkGmsh(ierr, "Mesh.MeshSizeExtendFromBoundary");
        gmshModelMeshSetSizeCallback(meshSizeAt, &grading, &ierr);
        checkGmsh(ierr, "gmshModelMeshSetSizeCallback");
    } else if (matrixSize > 0) {
        gmshOptionSetNumber("Mesh.MeshSizeMax", matrixSize, &ierr);
        checkGmsh(ierr, "Mesh.MeshSizeMax");
    }

    int meshError = 0;
    gmshModelMeshGenerate(3, &meshError);
    if (graded) {
        gmshModelMeshRemoveSizeCallback(&ierr);
        checkGmsh(ierr, "gmshModelMeshRemoveSizeCallback");
    }
    checkGmsh(meshError, "gmshModelMeshGenerate");
}

Mesh meshBox(const double *bounds, const mxArray *polygons, const MeshSizing &sizing)
{
    int ierr = 0;
    gmshOptionSetNumber("General.Terminal", 0, &ierr);
    checkGmsh(ierr, "gmshOptionSetNumber");
    gmshModelAdd("fracturedBox", &ierr);
    checkGmsh(ierr, "gmshModelAdd");
    const int box = gmshModelOccAddBox(bounds[0], bounds[1], bounds[2],
        bounds[3] - bounds[0], bounds[4] - bounds[1], bounds[5] - bounds[2], -1, &ierr);
    checkGmsh(ierr, "gmshModelOccAddBox");

    const size_t nFractures = mxGetNumberOfElements(polygons);
    // Preserve the input fracture identity when OCC splits a polygon.
    std::vector<std::vector<int>> fractureSurfaces(nFractures);
    FragmentMap fragments;
    if (nFractures) {
        std::vector<int> tools;
        for (size_t i = 0; i < nFractures; ++i) {
            tools.push_back(2);
            tools.push_back(addPolygon(mxGetCell(polygons, i)));
        }
        const int volume[] = {3, box};
        GmshBuffer<int> entities;
        gmshModelOccFragment(volume, 2, tools.data(), tools.size(),
            &entities.data, &entities.size, &fragments.entities.data, &fragments.sizes,
            &fragments.entities.size, -1, 1, 1, &ierr);
        checkGmsh(ierr, "gmshModelOccFragment");
        if (fragments.entities.size != nFractures + 1)
            throw std::runtime_error("Unexpected OCC fragment map size.");
        for (size_t i = 0; i < nFractures; ++i) {
            const int *entities = fragments.entities.data[i + 1];
            for (size_t j = 0; j < fragments.sizes[i + 1]; j += 2) {
                if (entities[j] == 2)
                    fractureSurfaces[i].push_back(entities[j + 1]);
            }
        }
    }

    gmshModelOccSynchronize(&ierr);
    checkGmsh(ierr, "gmshModelOccSynchronize");
    generateMesh(sizing, polygons);

    GmshBuffer<size_t> tags;
    GmshBuffer<double> coordinates, parametric;
    gmshModelMeshGetNodes(&tags.data, &tags.size, &coordinates.data, &coordinates.size,
                         &parametric.data, &parametric.size, -1, -1, 0, 0, &ierr);
    checkGmsh(ierr, "gmshModelMeshGetNodes");
    Mesh mesh;
    if (!tags.size || coordinates.size != 3 * tags.size)
        throw std::runtime_error("Gmsh returned invalid node coordinates.");
    mesh.coordinates.assign(coordinates.data, coordinates.data + coordinates.size);

    // Gmsh node tags need not be contiguous. Return MATLAB row indices.
    std::unordered_map<size_t, size_t> nodeRows;
    nodeRows.reserve(tags.size);
    for (size_t i = 0; i < tags.size; ++i)
        nodeRows.emplace(tags.data[i], i + 1);
    const int tet = gmshModelMeshGetElementType("Tetrahedron", 1, 0, &ierr);
    checkGmsh(ierr, "gmshModelMeshGetElementType");
    appendElements(tet, -1, nodeRows, mesh.tetrahedra);
    if (mesh.tetrahedra.empty())
        throw std::runtime_error("Gmsh returned no tetrahedra.");

    mesh.fractures.resize(nFractures);
    if (nFractures) {
        const int tri = gmshModelMeshGetElementType("Triangle", 1, 0, &ierr);
        checkGmsh(ierr, "gmshModelMeshGetElementType");
        for (size_t i = 0; i < nFractures; ++i) {
            for (int surface : fractureSurfaces[i])
                appendElements(tri, surface, nodeRows, mesh.fractures[i]);
            if (mesh.fractures[i].empty())
                throw std::runtime_error("Gmsh returned no triangles for fracture " +
                                         std::to_string(i + 1) + ".");
        }
    }
    return mesh;
}

Mesh runGmsh(const double *bounds, const mxArray *polygons, const MeshSizing &sizing)
{
    int ierr = 0;
    const int active = gmshIsInitialized(&ierr);
    checkGmsh(ierr, "gmshIsInitialized");
    if (active)
        throw std::runtime_error("Expected an inactive Gmsh session.");

    gmshInitialize(0, nullptr, 0, 0, &ierr);
    const bool initialized = ierr == 0;
    Mesh mesh;
    try {
        checkGmsh(ierr, "gmshInitialize");
        mesh = meshBox(bounds, polygons, sizing);
    } catch (...) {
        if (initialized || gmshIsInitialized(&ierr))
            gmshFinalize(&ierr);
        checkGmsh(ierr, "Gmsh session cleanup");
        throw;
    }
    gmshFinalize(&ierr);
    checkGmsh(ierr, "gmshFinalize");
    return mesh;
}

mxArray *toMatrix(const std::vector<double> &data, size_t width)
{
    const size_t rows = data.size() / width;
    mxArray *matrix = mxCreateDoubleMatrix(rows, width, mxREAL);
    double *out = mxGetPr(matrix);
    for (size_t row = 0; row < rows; ++row)
        for (size_t col = 0; col < width; ++col)
            out[col * rows + row] = data[row * width + col];
    return matrix;
}

bool isRealMatrix(const mxArray *array)
{
    return array && mxIsDouble(array) && !mxIsComplex(array) && !mxIsSparse(array) &&
           mxGetNumberOfDimensions(array) == 2;
}

} // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 5 || nlhs != 3 || !isRealMatrix(prhs[0]) ||
        mxGetM(prhs[0]) != 3 || mxGetN(prhs[0]) != 2 || !mxIsCell(prhs[1]))
        mexErrMsgIdAndTxt("GReS:GmshUsage",
            "Use [xyz, tets, triangles] = gmshMex(bounds, polygons, matrixSize, fractureSize, transitionFactor).");

    for (int arg = 2; arg < 5; ++arg) {
        if (!isRealMatrix(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1 ||
            !std::isfinite(mxGetScalar(prhs[arg])) || mxGetScalar(prhs[arg]) < 0)
            mexErrMsgIdAndTxt("GReS:GmshInput", "Expected finite, nonnegative sizing parameters.");
    }
    const MeshSizing sizing{mxGetScalar(prhs[2]), mxGetScalar(prhs[3]), mxGetScalar(prhs[4])};
    if (sizing.transition == 0)
        mexErrMsgIdAndTxt("GReS:GmshInput", "The transition factor must be positive.");

    const double *bounds = mxGetPr(prhs[0]);
    for (size_t axis = 0; axis < 3; ++axis) {
        if (!std::isfinite(bounds[axis]) || !std::isfinite(bounds[axis + 3]) ||
            bounds[axis] >= bounds[axis + 3])
            mexErrMsgIdAndTxt("GReS:GmshInput", "Invalid box bounds.");
    }

    const mxArray *polygons = prhs[1];
    for (size_t i = 0; i < mxGetNumberOfElements(polygons); ++i) {
        const mxArray *polygon = mxGetCell(polygons, i);
        if (!isRealMatrix(polygon) || mxGetM(polygon) < 3 || mxGetN(polygon) != 3)
            mexErrMsgIdAndTxt("GReS:GmshInput", "Expected an N-by-3 polygon array.");

        const double *xyz = mxGetPr(polygon);
        const size_t n = mxGetM(polygon);
        for (size_t j = 0; j < 3 * n; ++j) {
            const size_t axis = j / n;
            if (!std::isfinite(xyz[j]) || xyz[j] < bounds[axis] || xyz[j] > bounds[axis + 3])
                mexErrMsgIdAndTxt("GReS:GmshInput", "Fracture vertices must lie in the box.");
        }
    }

    try {
        const Mesh mesh = runGmsh(bounds, polygons, sizing);
        plhs[0] = toMatrix(mesh.coordinates, 3);
        plhs[1] = toMatrix(mesh.tetrahedra, 4);
        plhs[2] = mxCreateCellMatrix(mesh.fractures.size(), 1);
        for (size_t i = 0; i < mesh.fractures.size(); ++i)
            mxSetCell(plhs[2], i, toMatrix(mesh.fractures[i], 3));
    } catch (const std::exception &error) {
        mexErrMsgIdAndTxt("GReS:GmshAPI", "%s", error.what());
    }
}
