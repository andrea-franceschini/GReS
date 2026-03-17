//----------------------------------------------------------------------------------------
// readGMSHmesh_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// Preserves the original #ifdef MEX_FUNCTION dual-mode structure:
//   - Compiled with -DMEX_FUNCTION (default): MATLAB MEX gateway
//   - Compiled without it:                    standalone CLI executable
//
// MATLAB signature:
//   [coord, elems, regions] = readGMSHmesh(fileName)
//
// Build command (MEX):
//   mex -R2018a -DMEX_FUNCTION readGMSHmesh_wrap.cpp
//
// Build command (standalone):
//   g++ -std=c++17 -o readGMSHmesh_standalone readGMSHmesh_wrap.cpp
//
// ALL FIXES APPLIED (MEX side only — standalone main() unchanged)
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgIdAndTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() is only valid for arithmetic T.
//         Strings use factory.createCharArray() instead.
//
// [NEW-1] mxGetString() + fixed 1024-byte buffer replaced with
//         CharArray::toAscii() — no truncation risk, any path length.
//         The original truncation-detection check is also removed as
//         it is no longer needed.
//
// [NEW-2] mxCreateNumericMatrix + (double*) mxGetData() for coord
//         replaced with factory.createArray<double>({nNodes, 3}) and
//         flat column-major indexing — semantics identical to original.
//
// [NEW-3] mxCreateNumericMatrix + (int*) mxGetData() for elems
//         replaced with factory.createArray<int32_t>({nElems, cols}).
//         Typed accessor removes the void* cast.
//
// [NEW-4] mxCreateStructMatrix + mxSetField + mxCreateString for regions
//         replaced with factory.createStructArray() + StructArray field
//         assignment. Each scalar int field uses factory.createScalar<int32_t>();
//         the name string field uses factory.createCharArray().
//         This is the most significant structural change in this file.
//
// [NEW-5] NULL checks on mxCreate* return values removed — factory methods
//         in the C++ MEX API throw on allocation failure rather than
//         returning NULL, so explicit checks are redundant.
//
// [NEW-6] coord.resize(0) / elems.resize(0) / regions.resize(0) replaced
//         with .clear() which is the idiomatic C++ equivalent and more
//         clearly expresses the intent of releasing the data.
//
// [NEW-7] std::size_t used for all array sizes and loop bounds.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "readGMSHmesh.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

//========================================================================================
// MEX mode — C++ MEX API (R2018a+)
//========================================================================================

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>
#include <string>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // [NEW-1] Extract filename — no fixed buffer, no truncation risk
        // -----------------------------------------------------------------------
        const CharArray   filenameArr = inputs[0];
        const std::string fileName    = filenameArr.toAscii();

        // -----------------------------------------------------------------------
        // Call the GMSH reader
        // -----------------------------------------------------------------------
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
        // [NEW-2] Build coord output: nNodes x 3 double matrix (column-major)
        //
        //   col 0 (x): indices [0         .. nNodes-1  ]
        //   col 1 (y): indices [nNodes     .. 2*nNodes-1]
        //   col 2 (z): indices [2*nNodes   .. 3*nNodes-1]
        //
        //   Flat index for row i, col j  →  i + j * nNodes
        //   Identical to the original: ptrCoord[j*nNodes + i]
        // -----------------------------------------------------------------------
        const std::size_t nNodes = coord.size();

        TypedArray<double> mxCoord =
            factory.createArray<double>({nNodes, 3});

        for (std::size_t i = 0; i < nNodes; ++i) {
            mxCoord[i              ] = coord[i].x;   // col 0
            mxCoord[i +     nNodes ] = coord[i].y;   // col 1
            mxCoord[i + 2 * nNodes ] = coord[i].z;   // col 2
        }
        coord.clear();   // [NEW-6]

        // -----------------------------------------------------------------------
        // [NEW-3] Build elems output: nElems x (MAX_NUM_VERTICES+3) int32 matrix
        //
        //   Column layout (same as original):
        //     col 0          : element ID
        //     col 1          : tag
        //     col 2          : n (number of vertices in this element)
        //     col 3..2+MAX   : vertex indices; 0-padded for unused slots
        //
        //   Flat column-major index: i + j * nElems
        //   TypedArray is zero-initialised by the factory — no explicit fill needed
        //   for the padding zeros (factory.createArray zero-fills for numeric types).
        // -----------------------------------------------------------------------
        const std::size_t nElems  = elems.size();
        const std::size_t numCols = static_cast<std::size_t>(MAX_NUM_VERTICES) + 3;

        TypedArray<int32_t> mxElems =
            factory.createArray<int32_t>({nElems, numCols});

        for (std::size_t i = 0; i < nElems; ++i) {
            std::size_t j = 0;
            mxElems[i + j * nElems] = static_cast<int32_t>(elems[i].ID);
            j++;
            mxElems[i + j * nElems] = static_cast<int32_t>(elems[i].tag);
            j++;
            mxElems[i + j * nElems] = static_cast<int32_t>(elems[i].n);
            j++;
            // Filled vertices
            for (int k = 0; k < elems[i].n; ++k, ++j)
                mxElems[i + j * nElems] = static_cast<int32_t>(elems[i].v[k]);
            // Remaining slots already zero — factory zero-initialises numeric arrays
        }
        elems.clear();   // [NEW-6]

        // -----------------------------------------------------------------------
        // [NEW-4] Build regions output: nFields x 1 struct array
        //         Fields: "dim" (int32 scalar), "ID" (int32 scalar), "name" (char)
        //
        //   Original used mxCreateStructMatrix + mxSetField + mxCreateString.
        //   Modern equivalent: factory.createStructArray() + StructArray field
        //   assignment using factory-created scalars and char arrays.
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
        regions.clear();   // [NEW-6]

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(mxCoord);
        outputs[1] = std::move(mxElems);
        outputs[2] = std::move(mxRegions);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("MATLAB:mxImportGMSHmesh:invalidInputs",
                       "Only one input is required.");

        if (outputs.size() > 3)
            throwError("MATLAB:mxImportGMSHmesh:undefinedOutputs",
                       "Too many output arguments.");

        // Accept both 'char array' and "string" scalar from MATLAB caller
        const ArrayType t = inputs[0].getType();
        if (t != ArrayType::CHAR && t != ArrayType::MATLAB_STRING)
            throwError("MATLAB:mxImportGMSHmesh:inputNotString",
                       "Input must be a string or char array.");
    }

    // [FIX-E] createCharArray for strings; routes through MATLAB error()
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

//========================================================================================
// Standalone mode — unchanged from original
//========================================================================================

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
