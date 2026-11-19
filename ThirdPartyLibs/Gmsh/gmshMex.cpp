#include "mex.h"

extern "C" {
#include "gmshc.h"
}

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <cstring>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
#include <cstdint>

#ifdef __linux__
#include <signal.h>
#endif

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
        if (entities.data) {
            for (size_t i = 0; i < entities.size; ++i)
                gmshFree(entities.data[i]);
        }
        gmshFree(sizes);
    }
};

void checkGmsh(int code, const char *call)
{
    if (code != 0)
        throw std::runtime_error(std::string(call) + " failed (error " +
                                 std::to_string(code) + ").");
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

struct Call {
    int nout, nin;
    mxArray **out;
    const mxArray **in;
};

void require(bool condition, const char *message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void checkArity(const Call &call, int inputs, int outputs)
{
    require(call.nin == inputs && call.nout == outputs,
            "Incorrect number of inputs or outputs for this command.");
}

bool isRealDouble(const mxArray *array)
{
    return array && mxIsDouble(array) && !mxIsComplex(array) && !mxIsSparse(array);
}

template<class T> T number(const mxArray *array)
{
    require(isRealDouble(array) && mxGetNumberOfElements(array) == 1,
            "Expected a real double scalar.");
    const double value = mxGetScalar(array);
    require(std::isfinite(value), "Expected a finite scalar.");
    if (std::numeric_limits<T>::is_integer)
        require(std::floor(value) == value && value >= std::numeric_limits<T>::lowest() &&
                value <= std::numeric_limits<T>::max(), "Integer argument out of range.");
    return static_cast<T>(value);
}

std::string stringArgument(const mxArray *array)
{
    require(mxIsChar(array) && mxGetM(array) == 1, "Expected a character row vector.");
    char *buffer = mxArrayToString(array);
    require(buffer, "Could not convert character vector.");
    std::string value(buffer);
    mxFree(buffer);
    return value;
}

std::vector<int> intVector(const mxArray *array)
{
    require(mxIsInt32(array) && !mxIsComplex(array) &&
            (mxIsEmpty(array) || mxGetM(array) == 1 || mxGetN(array) == 1),
            "Expected an int32 vector.");
    if (mxIsEmpty(array))
        return {};
    const auto *data = static_cast<const int32_t *>(mxGetData(array));
    return {data, data + mxGetNumberOfElements(array)};
}

template<class T> mxArray *numericRow(const T *data, size_t size, mxClassID type)
{
    require(size == 0 || data, "Gmsh returned a missing output buffer.");
    mxArray *array = mxCreateNumericMatrix(1, size, type, mxREAL);
    if (size)
        std::memcpy(mxGetData(array), data, size * sizeof(T));
    return array;
}

// Instantiate a handler from each registered scalar Gmsh C signature.
// The final argument in these signatures is always the Gmsh error pointer.
template<class Signature, Signature Function> struct ScalarCommand;

template<class Result, class... Args, Result (*Function)(Args...)>
struct ScalarCommand<Result (*)(Args...), Function> {
    using Arguments = std::tuple<Args...>;

    template<size_t... I>
    static void invoke(const Call &call, const char *name, std::index_sequence<I...>)
    {
        int ierr = 0;
        const Result result = Function(number<typename std::tuple_element<I, Arguments>::type>(
                                           call.in[I])..., &ierr);
        checkGmsh(ierr, name);
        call.out[0] = mxCreateDoubleScalar(result);
    }

    static void run(const Call &call, const char *name)
    {
        constexpr size_t inputs = sizeof...(Args) - 1;
        checkArity(call, inputs, 1);
        invoke(call, name, std::make_index_sequence<inputs>{});
    }
};

#ifdef __linux__
// The Linux SDK resets these process-wide handlers on its first initialize.
// MATLAB's JVM relies on its SIGSEGV handler for normal execution.
class PreserveSignalHandlers {
public:
    PreserveSignalHandlers()
    {
        for (size_t i = 0; i < signals.size(); ++i)
            if (sigaction(signals[i], nullptr, &previous[i]) != 0)
                throw std::runtime_error("Could not save MATLAB signal handlers.");
    }

    ~PreserveSignalHandlers()
    {
        for (size_t i = 0; i < signals.size(); ++i)
            sigaction(signals[i], &previous[i], nullptr);
    }

private:
    const std::array<int, 11> signals{
        SIGHUP, SIGQUIT, SIGILL, SIGTRAP, SIGBUS, SIGFPE, SIGSEGV,
        SIGPIPE, SIGTERM, SIGURG, SIGSYS
    };
    std::array<struct sigaction, 11> previous{};
};
#endif

bool ownsSession = false;

void initialize(const Call &call)
{
    checkArity(call, 0, 0);
    int ierr = 0;
    const int active = gmshIsInitialized(&ierr);
    checkGmsh(ierr, "gmshIsInitialized");
    require(!active, "Expected an inactive Gmsh session.");
    {
#ifdef __linux__
        const PreserveSignalHandlers handlers;
#endif
        gmshInitialize(0, nullptr, 0, 0, &ierr);
    }
    if (ierr) {
        const int initializeError = ierr;
        if (gmshIsInitialized(&ierr))
            gmshFinalize(&ierr);
        checkGmsh(initializeError, "gmshInitialize");
    }
    ownsSession = true;
    mexLock();
}

void finalize(const Call &call)
{
    checkArity(call, 0, 0);
    int ierr = 0;
    gmshFinalize(&ierr);
    checkGmsh(ierr, "gmshFinalize");
    ownsSession = false;
    mexUnlock();
}

void addModel(const Call &call)
{
    checkArity(call, 1, 0);
    const auto name = stringArgument(call.in[0]);
    int ierr = 0;
    gmshModelAdd(name.c_str(), &ierr);
    checkGmsh(ierr, "gmshModelAdd");
}

void setNumber(const Call &call)
{
    checkArity(call, 2, 0);
    const auto name = stringArgument(call.in[0]);
    const double value = number<double>(call.in[1]);
    int ierr = 0;
    gmshOptionSetNumber(name.c_str(), value, &ierr);
    checkGmsh(ierr, "gmshOptionSetNumber");
}

void addOccLoop(const Call &call, bool surface)
{
    checkArity(call, 2, 1);
    const auto tags = intVector(call.in[0]);
    const int tag = number<int>(call.in[1]);
    const char *name = surface ? "gmshModelOccAddPlaneSurface" : "gmshModelOccAddCurveLoop";
    int ierr = 0;
    const int result = surface ? gmshModelOccAddPlaneSurface(tags.data(), tags.size(), tag, &ierr) :
                                 gmshModelOccAddCurveLoop(tags.data(), tags.size(), tag, &ierr);
    checkGmsh(ierr, name);
    call.out[0] = mxCreateDoubleScalar(result);
}

void fragment(const Call &call)
{
    checkArity(call, 5, 2);
    const auto objects = intVector(call.in[0]);
    const auto tools = intVector(call.in[1]);
    require(objects.size() % 2 == 0 && tools.size() % 2 == 0,
            "Expected alternating dimension and tag pairs.");
    const int tag = number<int>(call.in[2]);
    const int removeObject = number<int>(call.in[3]);
    const int removeTool = number<int>(call.in[4]);
    GmshBuffer<int> entities;
    FragmentMap map;
    int ierr = 0;
    gmshModelOccFragment(objects.data(), objects.size(), tools.data(), tools.size(),
                         &entities.data, &entities.size, &map.entities.data, &map.sizes,
                         &map.entities.size, tag, removeObject, removeTool, &ierr);
    checkGmsh(ierr, "gmshModelOccFragment");
    require(map.entities.size == 0 || (map.entities.data && map.sizes),
            "Gmsh returned an incomplete fragment map.");
    call.out[0] = numericRow(entities.data, entities.size, mxINT32_CLASS);
    call.out[1] = mxCreateCellMatrix(1, map.entities.size);
    for (size_t i = 0; i < map.entities.size; ++i)
        mxSetCell(call.out[1], i, numericRow(map.entities.data[i], map.sizes[i], mxINT32_CLASS));
}

void synchronize(const Call &call)
{
    checkArity(call, 0, 0);
    int ierr = 0;
    gmshModelOccSynchronize(&ierr);
    checkGmsh(ierr, "gmshModelOccSynchronize");
}

void generate(const Call &call)
{
    checkArity(call, 1, 0);
    const int dim = number<int>(call.in[0]);
    int ierr = 0;
    gmshModelMeshGenerate(dim, &ierr);
    checkGmsh(ierr, "gmshModelMeshGenerate");
}

void generateGraded(const Call &call)
{
    checkArity(call, 4, 0);
    const mxArray *polygons = call.in[0];
    require(mxIsCell(polygons), "Expected a cell array of polygons.");
    for (size_t i = 0; i < mxGetNumberOfElements(polygons); ++i) {
        const mxArray *polygon = mxGetCell(polygons, i);
        require(isRealDouble(polygon) && mxGetNumberOfDimensions(polygon) == 2 &&
                mxGetM(polygon) >= 3 && mxGetN(polygon) == 3,
                "Expected an N-by-3 real polygon array.");
    }
    const MeshSizing sizing{number<double>(call.in[1]), number<double>(call.in[2]),
                            number<double>(call.in[3])};
    require(sizing.matrix >= 0 && sizing.fracture >= 0 && sizing.transition > 0,
            "Invalid mesh sizing parameters.");
    generateMesh(sizing, polygons);
}

void getNodes(const Call &call)
{
    checkArity(call, 0, 2);
    GmshBuffer<size_t> tags;
    GmshBuffer<double> coordinates, parametric;
    int ierr = 0;
    gmshModelMeshGetNodes(&tags.data, &tags.size, &coordinates.data, &coordinates.size,
                         &parametric.data, &parametric.size, -1, -1, 0, 0, &ierr);
    checkGmsh(ierr, "gmshModelMeshGetNodes");
    require(coordinates.size == 3 * tags.size, "Invalid node coordinates returned by Gmsh.");
    call.out[0] = numericRow(tags.data, tags.size, mxUINT64_CLASS);
    call.out[1] = numericRow(coordinates.data, coordinates.size, mxDOUBLE_CLASS);
}

void getElementType(const Call &call)
{
    checkArity(call, 3, 1);
    const auto family = stringArgument(call.in[0]);
    const int order = number<int>(call.in[1]);
    const int serendip = number<int>(call.in[2]);
    int ierr = 0;
    const int type = gmshModelMeshGetElementType(family.c_str(), order, serendip, &ierr);
    checkGmsh(ierr, "gmshModelMeshGetElementType");
    call.out[0] = mxCreateDoubleScalar(type);
}

void getElementsByType(const Call &call)
{
    checkArity(call, 2, 2);
    const int type = number<int>(call.in[0]);
    const int tag = number<int>(call.in[1]);
    GmshBuffer<size_t> elements, nodes;
    int ierr = 0;
    gmshModelMeshGetElementsByType(type, &elements.data, &elements.size,
                                  &nodes.data, &nodes.size, tag, 0, 1, &ierr);
    checkGmsh(ierr, "gmshModelMeshGetElementsByType");
    call.out[0] = numericRow(elements.data, elements.size, mxUINT64_CLASS);
    call.out[1] = numericRow(nodes.data, nodes.size, mxUINT64_CLASS);
}

} // namespace

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    try {
        require(nrhs >= 1, "Specify a Gmsh command.");
        const auto command = stringArgument(prhs[0]);
        const Call call{nlhs, nrhs - 1, plhs, prhs + 1};
        if (command == "initialize") {
            initialize(call);
            return;
        }
        require(ownsSession, "Initialize Gmsh before calling this command.");

        if (command == "finalize")
            finalize(call);
        else if (command == "model.add")
            addModel(call);
        else if (command == "option.setNumber")
            setNumber(call);
        else if (command == "model.occ.addBox")
            ScalarCommand<decltype(&gmshModelOccAddBox), gmshModelOccAddBox>::run(call, "gmshModelOccAddBox");
        else if (command == "model.occ.addPoint")
            ScalarCommand<decltype(&gmshModelOccAddPoint), gmshModelOccAddPoint>::run(call, "gmshModelOccAddPoint");
        else if (command == "model.occ.addLine")
            ScalarCommand<decltype(&gmshModelOccAddLine), gmshModelOccAddLine>::run(call, "gmshModelOccAddLine");
        else if (command == "model.occ.addCurveLoop")
            addOccLoop(call, false);
        else if (command == "model.occ.addPlaneSurface")
            addOccLoop(call, true);
        else if (command == "model.occ.fragment")
            fragment(call);
        else if (command == "model.occ.synchronize")
            synchronize(call);
        else if (command == "model.mesh.generate")
            generate(call);
        else if (command == "model.mesh.generateGraded")
            generateGraded(call);
        else if (command == "model.mesh.getNodes")
            getNodes(call);
        else if (command == "model.mesh.getElementType")
            getElementType(call);
        else if (command == "model.mesh.getElementsByType")
            getElementsByType(call);
        else
            throw std::runtime_error("Unsupported Gmsh command: " + command);
    } catch (const std::exception &error) {
        mexErrMsgIdAndTxt("GReS:GmshAPI", "%s", error.what());
    }
}
