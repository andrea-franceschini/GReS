//----------------------------------------------------------------------------------------
// readGMSHmesh_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// Preserves the original #ifdef MEX_FUNCTION dual-mode structure.
//
// MATLAB signature:
//   [coord, elems, regions] = readGMSHmesh(fileName)
//
// Build command (MEX):
//   mex -R2018a -DMEX_FUNCTION readGMSHmesh_wrap.cpp
//
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-G] TypedArray flat-index subscripting causes runtime error
//         "Not enough indices provided" on 2D arrays.
//
//         mxCoord {nNodes, 3}:
//           mxCoord[i + j*nNodes] was flat-indexing a 2D TypedArray.
//           Fixed: fill std::vector<double> coordBuf column-major,
//           then std::copy into mxCoord via iterators.
//
//         mxElems {nElems, numCols}:
//           mxElems[i + j*nElems] was flat-indexing a 2D TypedArray.
//           Fixed: fill std::vector<int32_t> elemsBuf column-major,
//           then std::copy into mxElems via iterators.
//
// All previously documented fixes (FIX-A/C/E, NEW-1..6) are retained.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "readGMSHmesh.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <string>
#include <algorithm> // std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        const CharArray   filenameArr = inputs[0];
        const std::string fileName    = filenameArr.toAscii();

        std::vector<Region>  regions;
        std::vector<Point>   coord;
        std::vector<Element> elems;

        int status = readGMSHmesh(fileName, coord, elems, regions);

        if (status == 1)
            throwError("MATLAB:mxImportGMSHmesh",
                       "Input file does not exist.");
        else if (status == 2)
            throwError("MATLAB:mxImportGMSHmesh",
                       "Wrong file format.");

        // -----------------------------------------------------------------------
        // [FIX-G] coord output: nNodes x 3 double
        //   Fill flat column-major buffer, then std::copy into TypedArray.
        // -----------------------------------------------------------------------
        const std::size_t nNodes = coord.size();

        std::vector<double> coordBuf(nNodes * 3);
        for (std::size_t i = 0; i < nNodes; ++i) {
            coordBuf[i              ] = coord[i].x;
            coordBuf[i +     nNodes ] = coord[i].y;
            coordBuf[i + 2 * nNodes ] = coord[i].z;
        }
        coord.clear();

        TypedArray<double> mxCoord = factory.createArray<double>({nNodes, 3});
        std::copy(coordBuf.begin(), coordBuf.end(), mxCoord.begin());

        // -----------------------------------------------------------------------
        // [FIX-G] elems output: nElems x (MAX_NUM_VERTICES+3) int32
        //   Fill flat column-major buffer (zero-init for padding),
        //   then std::copy into TypedArray.
        // -----------------------------------------------------------------------
        const std::size_t nElems  = elems.size();
        const std::size_t numCols = static_cast<std::size_t>(MAX_NUM_VERTICES) + 3;

        std::vector<int32_t> elemsBuf(nElems * numCols, 0);  // zero-init = padding

        for (std::size_t i = 0; i < nElems; ++i) {
            std::size_t j = 0;
            elemsBuf[i + j * nElems] = static_cast<int32_t>(elems[i].ID);   ++j;
            elemsBuf[i + j * nElems] = static_cast<int32_t>(elems[i].tag);  ++j;
            elemsBuf[i + j * nElems] = static_cast<int32_t>(elems[i].n);    ++j;
            for (int k = 0; k < elems[i].n; ++k, ++j)
                elemsBuf[i + j * nElems] = static_cast<int32_t>(elems[i].v[k]);
            // remaining slots already 0 from zero-init
        }
        elems.clear();

        TypedArray<int32_t> mxElems =
            factory.createArray<int32_t>({nElems, numCols});
        std::copy(elemsBuf.begin(), elemsBuf.end(), mxElems.begin());

        // -----------------------------------------------------------------------
        // regions output: nFields x 1 struct array (no flat indexing — safe)
        // -----------------------------------------------------------------------
        const std::size_t nFields = regions.size();

        StructArray mxRegions =
            factory.createStructArray({nFields, 1}, {"dim", "ID", "name"});

        for (std::size_t i = 0; i < nFields; ++i) {
            mxRegions[i]["dim"]  = factory.createScalar<int32_t>(
                                       static_cast<int32_t>(regions[i].dim));
            mxRegions[i]["ID"]   = factory.createScalar<int32_t>(
                                       static_cast<int32_t>(regions[i].ID));
            mxRegions[i]["name"] = factory.createCharArray(regions[i].name);
        }
        regions.clear();

        outputs[0] = std::move(mxCoord);
        outputs[1] = std::move(mxElems);
        outputs[2] = std::move(mxRegions);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("MATLAB:mxImportGMSHmesh:invalidInputs",
                       "Only one input is required.");

        if (outputs.size() > 3)
            throwError("MATLAB:mxImportGMSHmesh:undefinedOutputs",
                       "Too many output arguments.");

        const ArrayType t = inputs[0].getType();
        if (t != ArrayType::CHAR && t != ArrayType::MATLAB_STRING)
            throwError("MATLAB:mxImportGMSHmesh:inputNotString",
                       "Input must be a string or char array.");
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

#else // MEX_FUNCTION

#include <iostream>
#include <vector>
#include <string>

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cout << "Only one input is required." << std::endl;
        return 1;
    }

    std::string const fileName(argv[1]);
    std::vector<Point>   coord;
    std::vector<Element> elems;
    std::vector<Region>  regions;
    readGMSHmesh(fileName, coord, elems, regions);

    for (int i = 0; i < (int)coord.size(); ++i)
        std::cout << coord[i].x << " "
                  << coord[i].y << " "
                  << coord[i].z << std::endl;

    for (int i = 0; i < (int)elems.size(); ++i) {
        std::cout << elems[i].ID << " " << elems[i].tag << " ";
        for (int j = 0; j < elems[i].n; ++j)
            std::cout << elems[i].v[j] << " ";
        std::cout << std::endl;
    }

    for (int i = 0; i < (int)regions.size(); ++i)
        std::cout << regions[i].dim << " "
                  << regions[i].ID  << " "
                  << regions[i].name << " " << std::endl;

    std::cout << "END" << std::endl;
    return 0;
}

#endif // MEX_FUNCTION

//----------------------------------------------------------------------------------------
