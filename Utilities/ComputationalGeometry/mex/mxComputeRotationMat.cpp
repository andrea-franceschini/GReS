//----------------------------------------------------------------------------------------
// computeRot_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   R = computeRot(n)
//     n : 3-element real double vector (normal)
//     R : 3 x 3 real double rotation matrix (columns: n, m1, m2)
//
// Build command:
//   mex -R2018a computeRot_wrap.cpp
//
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-G] TypedArray flat-index subscripting causes runtime error
//         "Not enough indices provided" on 2D arrays.
//         The previous version wrote R_out[k] in a loop over a {3,3}
//         TypedArray — single flat index on a 2D array is illegal.
//         Fixed by filling a std::vector<double>(9) column-major and
//         then std::copy-ing into R_out via iterators.
//
// All previously documented fixes (FIX-A/B/C/E, NEW-1..4) are retained.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm> // std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

//========================================================================================
// Pure math helpers — no MATLAB API dependency
//========================================================================================

inline std::array<double,3> cross(const std::array<double,3>& a,
                                  const std::array<double,3>& b)
{
    return { a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] };
}

inline double dot(const std::array<double,3>& a, const std::array<double,3>& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double norm(const std::array<double,3>& a)
{
    return std::sqrt(dot(a, a));
}

inline std::array<double,3> normalize(const std::array<double,3>& a)
{
    const double n = norm(a);
    if (n < 1e-16)
        throw std::runtime_error("Zero-length vector cannot be normalized.");
    return { a[0]/n, a[1]/n, a[2]/n };
}

inline double det3x3(const std::array<std::array<double,3>,3>& R)
{
    const auto& a = R[0]; const auto& b = R[1]; const auto& c = R[2];
    return a[0]*(b[1]*c[2] - b[2]*c[1])
         - a[1]*(b[0]*c[2] - b[2]*c[0])
         + a[2]*(b[0]*c[1] - b[1]*c[0]);
}

//========================================================================================
// MexFunction
//========================================================================================
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        // Read n
        const TypedArray<double> nArr = inputs[0];
        std::array<double,3> n = { nArr[0], nArr[1], nArr[2] };

        try {
            n = normalize(n);

            // Pick a vector not parallel to n
            std::array<double,3> tmp;
            if (std::fabs(n[0]) < 0.9) tmp = {1.0, 0.0, 0.0};
            else                        tmp = {0.0, 1.0, 0.0};

            // First tangent
            const double ndot = dot(tmp, n);
            std::array<double,3> m1 = { tmp[0] - ndot*n[0],
                                        tmp[1] - ndot*n[1],
                                        tmp[2] - ndot*n[2] };
            m1 = normalize(m1);

            // Second tangent
            std::array<double,3> m2 = cross(n, m1);
            m2 = normalize(m2);

            // Assemble columns: n, m1, m2
            std::array<std::array<double,3>,3> Rcols = { n, m1, m2 };

            // Enforce right-handed system
            double detR = det3x3(Rcols);
            if (detR < 0.0) {
                m1[0] = -m1[0]; m1[1] = -m1[1]; m1[2] = -m1[2];
                Rcols = { n, m1, m2 };
                detR  = det3x3(Rcols);
            }

            if (std::fabs(detR - 1.0) > 1e-12)
                throwError("computeRot:notOrthogonal",
                           "Rotation matrix is not orthogonal to machine precision.");

            // -------------------------------------------------------------------
            // [FIX-G] Build 3x3 output via flat std::vector (column-major),
            //         then std::copy into TypedArray via iterators.
            //         Direct R_out[k] subscripting on a {3,3} TypedArray
            //         raises "Not enough indices provided" at runtime.
            // -------------------------------------------------------------------
            std::vector<double> Rbuf(9);
            for (int j = 0; j < 3; ++j)
                for (int i = 0; i < 3; ++i)
                    Rbuf[i + j*3] = Rcols[j][i];   // column-major

            TypedArray<double> R_out = factory.createArray<double>({3, 3});
            std::copy(Rbuf.begin(), Rbuf.end(), R_out.begin());

            outputs[0] = std::move(R_out);

        } catch (const std::runtime_error& e) {
            throwError("computeRot:mathError", e.what());
        }
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("computeRot:inputError",
                       "Usage: R = computeRot(n).");

        if (outputs.size() > 1)
            throwError("computeRot:outputError",
                       "One output required.");

        if (inputs[0].getType() != ArrayType::DOUBLE)
            throwError("computeRot:inputError",
                       "Input n must be a 3-element real vector.");

        if (inputs[0].getNumberOfElements() != 3)
            throwError("computeRot:inputError",
                       "Input n must be a 3-element real vector.");
    }

    void throwError(const std::string& id, const std::string& msg)
    {
        getEngine()->feval(u"error", 0,
            std::vector<Array>{
                factory.createCharArray(id),
                factory.createCharArray(msg)
            });
    }
};

//----------------------------------------------------------------------------------------
