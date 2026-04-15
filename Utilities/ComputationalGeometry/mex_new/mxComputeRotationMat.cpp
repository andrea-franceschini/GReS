#include "mex.hpp"
#include "mexAdapter.hpp"
#include "src/MexHelper.hpp"

#include <vector>

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
        requireArgs(this, factory, inputs.size() == 1,
                    "mxComputeRotationMat:input", "Usage: R = mxComputeRotationMat(normal).");
        requireArgs(this, factory, outputs.size() <= 1,
                    "mxComputeRotationMat:output", "One output only.");
        requireArgs(this, factory,
                    inputs[0].getType() == ArrayType::DOUBLE && numel(inputs[0]) == 3,
                    "mxComputeRotationMat:input", "normal must be a real 3-vector.");

        std::vector<double> nbuf = typedArrayToVector(inputs[0]);

        std::vector<double> R(9, 0.0);
        try {
            polygeom::rotationFromNormal(nbuf.data(), R.data());
        } catch (const std::exception& e) { throwErr("mxComputeRotationMat:error", e.what()); }

        outputs[0] = vectorToArray(factory, R, 3, 3);
    }
};