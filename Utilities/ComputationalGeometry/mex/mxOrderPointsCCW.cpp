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

    TypedArray<double> localCall(const Array& pointsArr, const Array* normalArr)
    {
        auto dims = pointsArr.getDimensions();
        std::size_t n = dims[0], dim = dims[1];
        requireArgs(this, factory, (dim == 2 || dim == 3) && n >= 3,
                    "mxOrderPointsCCW:input", "points must be N x 2 or N x 3 with N >= 3.");

        std::vector<double> P = typedArrayToVector(pointsArr);
        std::vector<polygeom::Index> perm;

        try {
            if (dim == 2) {
                polygeom::orderCCW2D(P.data(), n, perm);
            } else {
                std::vector<double> normalBuf;
                const double* userNormal = nullptr;
                if (normalArr) {
                    requireArgs(this, factory,
                                (*normalArr).getType() == ArrayType::DOUBLE && numel(*normalArr) == 3,
                                "mxOrderPointsCCW:input", "normal must be a real 3-vector.");
                    normalBuf = typedArrayToVector(*normalArr);
                    userNormal = normalBuf.data();
                }
                polygeom::orderCCW3D(P.data(), n, userNormal, perm, nullptr);
            }
        } catch (const std::exception& e) { throwErr("mxOrderPointsCCW:error", e.what()); }

        std::vector<double> idx(n);
        for (std::size_t i = 0; i < n; ++i) idx[i] = static_cast<double>(perm[i] + 1);
        return vectorToArray(factory, idx, n, 1);
    }

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        requireArgs(this, factory, outputs.size() <= 2,
                    "mxOrderPointsCCW:output", "At most two outputs: [idx] or [Pccw, perm].");
        requireArgs(this, factory, inputs.size() >= 1 && inputs.size() <= 3,
                    "mxOrderPointsCCW:input",
                    "Usage: idx = mxOrderPointsCCW(P [,normal]) or [Pccw,perm] = mxOrderPointsCCW(Pflat,nVert [,normals]).");

        bool batchCall = false;
        if (inputs.size() >= 2) {
            bool secondIsVector = (getRows(inputs[1]) == 1 || getCols(inputs[1]) == 1);
            bool firstLooksMatrix = (getCols(inputs[0]) == 2 || getCols(inputs[0]) == 3);
            bool secondHasMany = (numel(inputs[1]) > 1);
            batchCall = firstLooksMatrix && secondIsVector && secondHasMany;
        }

        if (!batchCall) {
            requireArgs(this, factory, outputs.size() <= 1,
                        "mxOrderPointsCCW:output", "Local call returns only idx.");
            const Array* normalArr = (inputs.size() >= 2) ? &inputs[1] : nullptr;
            outputs[0] = localCall(inputs[0], normalArr);
            return;
        }

        std::vector<double> Pbuf, nVertBuf;
        polygeom::BatchInput in;
        try {
            in = parseBatchInput(this, factory, inputs[0], inputs[1], Pbuf, nVertBuf);
        } catch (const std::exception& e) { throwErr("mxOrderPointsCCW:input", e.what()); }

        std::vector<double> normalsBuf;
        const double* normalsOrNull = nullptr;
        if (inputs.size() == 3) {
            requireArgs(this, factory, in.dim == 3,
                        "mxOrderPointsCCW:input",
                        "Batch normals can only be supplied for 3D polygons.");
            requireArgs(this, factory,
                        inputs[2].getType() == ArrayType::DOUBLE &&
                        getRows(inputs[2]) == in.nPoly && getCols(inputs[2]) == 3,
                        "mxOrderPointsCCW:input", "normals must have size nPoly x 3.");
            normalsBuf = typedArrayToVector(inputs[2]);
            normalsOrNull = normalsBuf.data();
        }

        std::vector<double> Pccw, perm;
        try {
            polygeom::orderPointsBatch(in, normalsOrNull, Pccw, perm);
        } catch (const std::exception& e) { throwErr("mxOrderPointsCCW:error", e.what()); }

        outputs[0] = vectorToArray(factory, Pccw, in.nPts, in.dim);
        if (outputs.size() >= 2)
            outputs[1] = vectorToArray(factory, perm, in.nPts, 1);
    }
};