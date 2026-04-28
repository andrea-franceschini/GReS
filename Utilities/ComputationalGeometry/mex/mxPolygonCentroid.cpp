#include "mex.hpp"
#include "mexAdapter.hpp"
#include "src/MexHelper.hpp"

#include <vector>
#include <algorithm>

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace mexhelper;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

    void throwErr(const std::string& id, const std::string& msg) {
        throwError(this, factory, id, msg);
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        requireArgs(this, factory, outputs.size() <= 1,
                    "mxPolygonCentroid:output", "One output only.");

        if (inputs.size() == 1) {
            validateDoubleMatrix(this, factory, inputs[0],
                                 "mxPolygonCentroid:input", "points");
            auto dims = inputs[0].getDimensions();
            std::size_t nPts = dims[0], dim = dims[1];
            requireArgs(this, factory, (dim == 2 || dim == 3) && nPts >= 3,
                        "mxPolygonCentroid:input", "points must be N x 2 or N x 3 with N >= 3.");

            std::vector<double> P = typedArrayToVector(inputs[0]);

            double area, c[3] = {0,0,0}, normal[3];
            try {
                if (dim == 2) polygeom::areaCentroidNormal2D(P.data(), nPts, area, c);
                else          polygeom::areaCentroidNormal3D(P.data(), nPts, nullptr, area, c, normal);
            } catch (const std::exception& e) { throwErr("mxPolygonCentroid:error", e.what()); }

            std::vector<double> cvec(c, c + dim);
            outputs[0] = vectorToArray(factory, cvec, 1, dim);
            return;
        }

        requireArgs(this, factory, inputs.size() == 2,
                    "mxPolygonCentroid:input",
                    "Usage: C = mxPolygonCentroid(points) or C = mxPolygonCentroid(Pflat,nVert).");

        std::vector<double> Pbuf, nVertBuf;
        polygeom::BatchInput in;
        try {
            in = parseBatchInput(this, factory, inputs[0], inputs[1], Pbuf, nVertBuf);
        } catch (const std::exception& e) { throwErr("mxPolygonCentroid:input", e.what()); }

        std::vector<double> centroid;
        try {
            polygeom::polygonCentroidBatch(in, nullptr, centroid);
        } catch (const std::exception& e) { throwErr("mxPolygonCentroid:error", e.what()); }

        outputs[0] = vectorToArray(factory, centroid, in.nPoly, in.dim);
    }
};