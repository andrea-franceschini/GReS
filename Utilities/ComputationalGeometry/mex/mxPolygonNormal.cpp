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
                    "mxPolygonNormal:output", "One output only.");

        if (inputs.size() == 1) {
            validateDoubleMatrix(this, factory, inputs[0],
                                 "mxPolygonNormal:input", "points");
            auto dims = inputs[0].getDimensions();
            std::size_t nPts = dims[0], dim = dims[1];
            requireArgs(this, factory, dim == 3 && nPts >= 3,
                        "mxPolygonNormal:input", "points must be N x 3 with N >= 3.");

            std::vector<double> P = typedArrayToVector(inputs[0]);

            double area, c[3], normal[3];
            try {
                polygeom::areaCentroidNormal3D(P.data(), nPts, nullptr, area, c, normal);
            } catch (const std::exception& e) { throwErr("mxPolygonNormal:error", e.what()); }

            std::vector<double> nvec(normal, normal + 3);
            outputs[0] = vectorToArray(factory, nvec, 1, 3);
            return;
        }

        requireArgs(this, factory, inputs.size() == 2,
                    "mxPolygonNormal:input",
                    "Usage: N = mxPolygonNormal(points3D) or N = mxPolygonNormal(Pflat3D,nVert).");

        std::vector<double> Pbuf, nVertBuf;
        polygeom::BatchInput in;
        try {
            in = parseBatchInput(this, factory, inputs[0], inputs[1], Pbuf, nVertBuf);
        } catch (const std::exception& e) { throwErr("mxPolygonNormal:input", e.what()); }

        std::vector<double> normal;
        try {
            polygeom::polygonNormalBatch(in, nullptr, normal);
        } catch (const std::exception& e) { throwErr("mxPolygonNormal:error", e.what()); }

        outputs[0] = vectorToArray(factory, normal, in.nPoly, 3);
    }
};