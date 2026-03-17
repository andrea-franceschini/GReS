//----------------------------------------------------------------------------------------
// VTUWriter_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   VTUWriter_wrap(fileName, time, coord, cells, cellVTKType,
//                  cellNumVerts [, fieldStruct1, fieldStruct2, ...])
//
// Build command:
//   mex -R2018a -DMEX_FUNCTION VTUWriter_wrap.cpp
//
// ALL FIXES APPLIED
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
// [NEW-1] mxGetString() + fixed char[] buffer replaced with
//         CharArray::toAscii() — no truncation risk, any path length.
//         Applied to both the fileName and all struct fieldName strings.
//
// [NEW-2] mxGetPr() on coord, cells, cellVTKType, cellNumVerts replaced
//         with TypedArray<double> views. Raw data is accessed via flat
//         iterator index preserving original column-major semantics.
//
// [NEW-3] malloc/free for cellNumVertsInt, cellVTKTypeInt, cellsInt
//         replaced with std::vector<int> — RAII, exception-safe.
//
// [NEW-4] Struct field processing (the most complex part):
//         mxGetClassID / mxGetNumberOfFields / mxGetField / mxIsChar /
//         mxGetData on struct inputs replaced with:
//           - inputs[i].getType() == ArrayType::STRUCT
//           - StructArray::getNumberOfFields()
//           - StructArray::getFieldNames() to verify "name" and "data" exist
//           - StructArray[i]["name"] → CharArray → toAscii()
//           - StructArray[i]["data"] → TypedArray<double>
//           - TypedArray::getDimensions() replaces mxGetM / mxGetN
//
// [NEW-5] mxGetM / mxGetN replaced with getDimensions()[0] / [1]
//         throughout for coord, cells, and field data arrays.
//
// [NEW-6] NULL checks on mxGetField replaced with field-name existence
//         checks via getFieldNames() before accessing field values.
//
// [NEW-7] std::size_t / int used consistently; no implicit narrowing casts.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "VTUWriter.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <vector>
#include <string>
#include <algorithm>  // std::find

using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// Helper: check whether a field name exists in a StructArray
//----------------------------------------------------------------------------------------
static bool hasField(const StructArray& sa, const std::string& name)
{
    for (const auto& fn : sa.getFieldNames())
        if (std::string(fn) == name) return true;
    return false;
}

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // [NEW-1] Extract filename
        // -----------------------------------------------------------------------
        const CharArray   fileNameArr = inputs[0];
        const std::string fileName    = fileNameArr.toAscii();

        // -----------------------------------------------------------------------
        // Read scalar: time
        // -----------------------------------------------------------------------
        const double time = static_cast<double>(
                                TypedArray<double>(inputs[1])[0]);

        // -----------------------------------------------------------------------
        // [NEW-2][NEW-5] Read coord — nPoints x dim double matrix
        // -----------------------------------------------------------------------
        const TypedArray<double> mxCoord = inputs[2];
        const auto coordDims = mxCoord.getDimensions();
        const int nPoints    = static_cast<int>(coordDims[0]);
        const int dim        = static_cast<int>(coordDims[1]);

        // Copy into contiguous vector so we can pass a raw pointer to VTUWriter
        std::vector<double> coordVec(mxCoord.begin(), mxCoord.end());

        // -----------------------------------------------------------------------
        // [NEW-2][NEW-5] Read cells, cellVTKType, cellNumVerts
        // -----------------------------------------------------------------------
        const TypedArray<double> mxCells        = inputs[3];
        const TypedArray<double> mxCellVTKType  = inputs[4];
        const TypedArray<double> mxCellNumVerts = inputs[5];

        const auto cellsDims       = mxCells.getDimensions();
        const int  nCells          = static_cast<int>(cellsDims[0]);
        const int  maxCellNumVerts = static_cast<int>(cellsDims[1]);

        // -----------------------------------------------------------------------
        // Count non-zero connections (same logic as original)
        // -----------------------------------------------------------------------
        int nConnections = 0;
        {
            std::size_t total = static_cast<std::size_t>(nCells) *
                                static_cast<std::size_t>(maxCellNumVerts);
            auto it = mxCells.begin();
            for (std::size_t idx = 0; idx < total; ++idx, ++it)
                if (static_cast<int>(*it) != 0) ++nConnections;
        }

        // -----------------------------------------------------------------------
        // [NEW-3] Build int conversion arrays as std::vector (replaces malloc)
        // -----------------------------------------------------------------------
        std::vector<int> cellNumVertsInt(static_cast<std::size_t>(nCells));
        std::vector<int> cellVTKTypeInt (static_cast<std::size_t>(nCells));

        for (int i = 0; i < nCells; ++i) {
            cellNumVertsInt[static_cast<std::size_t>(i)] =
                static_cast<int>(mxCellNumVerts[static_cast<std::size_t>(i)]);
            cellVTKTypeInt [static_cast<std::size_t>(i)] =
                static_cast<int>(mxCellVTKType [static_cast<std::size_t>(i)]);
        }

        // Build flat connectivity array, converting MATLAB 1-based → 0-based
        std::vector<int> cellsInt(static_cast<std::size_t>(nConnections));
        {
            int k = 0;
            for (int i = 0; i < nCells; ++i) {
                for (int j = 0; j < cellNumVertsInt[static_cast<std::size_t>(i)]; ++j) {
                    // Column-major flat index: i + nCells*j
                    const std::size_t idx =
                        static_cast<std::size_t>(i) +
                        static_cast<std::size_t>(nCells) *
                        static_cast<std::size_t>(j);
                    cellsInt[static_cast<std::size_t>(k++)] =
                        static_cast<int>(mxCells[idx]) - 1;  // 1-based → 0-based
                }
            }
        }

        // -----------------------------------------------------------------------
        // Validate cell types (unchanged logic)
        // -----------------------------------------------------------------------
        for (int i = 0; i < nCells; ++i) {
            const int t = cellVTKTypeInt[static_cast<std::size_t>(i)];
            if (t != VTK_TRIANGLE    && t != VTK_QUAD      &&
                t != VTK_TETRA       && t != VTK_HEXAHEDRON &&
                t != VTK_HEXAHEDRON2 && t != VTK_QUAD2      &&
                t != VTK_POLYGON)
                throwError("MATLAB:mxVTKWriter:cellType",
                           "Element type not yet supported.");
        }

        // -----------------------------------------------------------------------
        // Set up VTUWriter and write mesh
        // -----------------------------------------------------------------------
        VTUWriter vtuFile;
        vtuFile.set_mesh(dim, nPoints, coordVec.data(),
                         nCells, nConnections,
                         cellsInt.data(),
                         cellVTKTypeInt.data(),
                         cellNumVertsInt.data());
        vtuFile.set_time(time);

        // -----------------------------------------------------------------------
        // [NEW-4] Process optional field structs (inputs 6, 7, ...)
        //
        //   Original:
        //     mxGetClassID / mxGetNumberOfFields / mxGetField /
        //     mxIsChar / mxGetData / mxGetM / mxGetN
        //
        //   Modern:
        //     getType() == ArrayType::STRUCT
        //     StructArray::getNumberOfFields()
        //     hasField() + StructArray[i]["name"] / ["data"]
        //     TypedArray<double>::getDimensions()
        // -----------------------------------------------------------------------
        const std::size_t nFixed = 6;
        for (std::size_t iStruct = 0; iStruct < inputs.size() - nFixed; ++iStruct) {

            Array inputArg = inputs[nFixed + iStruct];

            if (inputArg.getNumberOfElements() == 0) continue;

            // [NEW-4] Type check
            if (inputArg.getType() != ArrayType::STRUCT)
                throwError("MATLAB:mxVTKWriter:struct",
                           "Scalar and vector data must be saved as struct.");

            StructArray sa = inputs[nFixed + iStruct];

            // [NEW-4] Field count check
            if (sa.getNumberOfFields() != 2)
                throwError("MATLAB:mxVTKWriter:struct",
                           "Struct for scalar and vector data must have exactly 2 fields.");

            // [NEW-4] Verify required field names exist before accessing
            if (!hasField(sa, "name"))
                throwError("MATLAB:mxVTKWriter:invalidField",
                           "name field must exist.");
            if (!hasField(sa, "data"))
                throwError("MATLAB:mxVTKWriter:invalidField",
                           "data field must exist.");

            const std::size_t nStructElems = sa.getNumberOfElements();

            for (std::size_t i = 0; i < nStructElems; ++i) {

                // [NEW-1][NEW-4] Extract field name string
                const Array nameRaw = sa[i]["name"];
                if (nameRaw.getType() != ArrayType::CHAR &&
                    nameRaw.getType() != ArrayType::MATLAB_STRING)
                    throwError("MATLAB:mxVTKWriter:invalidField",
                               "name field must be a string.");

                const CharArray   nameArr  = sa[i]["name"];
                const std::string fieldName = nameArr.toAscii();

                // [NEW-4] Extract data field
                const Array dataRaw = sa[i]["data"];
                if (dataRaw.getType() != ArrayType::DOUBLE)
                    throwError("MATLAB:mxVTKWriter:invalidField",
                               "data field must be a double array.");

                const TypedArray<double> dataArr = sa[i]["data"];
                const auto dataDims = dataArr.getDimensions();    // [NEW-5]
                const int  dataRows = static_cast<int>(dataDims[0]);
                const int  dataCols = static_cast<int>(dataDims[1]);

                if (dataRows != nPoints && dataRows != nCells)
                    throwError("MATLAB:mxVTKWriter:invalidField",
                               "Size inconsistency in data field.");

                // Copy data into contiguous vector for VTUWriter raw pointer call
                std::vector<double> dataVec(dataArr.begin(), dataArr.end());

                if (dataRows == nPoints)
                    vtuFile.add_field<double>(fieldName,
                                              dataVec.data(), dataRows, dataCols);
                else
                    vtuFile.add_cell_field<double>(fieldName,
                                                   dataVec.data(), dataRows, dataCols);
            }
        }

        // Write file — vectors (cellsInt etc.) self-destruct after this [NEW-3]
        vtuFile.write_mesh(fileName.c_str());
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() < 7)
            throwError("MATLAB:mxVTKWriter:invalidInputs",
                       "At least 7 input arguments are required.");

        if (outputs.size() > 0)
            throwError("MATLAB:mxVTKWriter:invalidOutputs",
                       "Too many output arguments.");

        // fileName must be a string
        const ArrayType t = inputs[0].getType();
        if (t != ArrayType::CHAR && t != ArrayType::MATLAB_STRING)
            throwError("MATLAB:mxVTKWriter:inputNotString",
                       "First input (fileName) must be a string or char array.");

        // time must be a scalar double
        if (inputs[1].getType() != ArrayType::DOUBLE ||
            inputs[1].getNumberOfElements() != 1)
            throwError("MATLAB:mxVTKWriter:invalidInputs",
                       "Second input (time) must be a real double scalar.");

        // coord, cells, cellVTKType, cellNumVerts must be double arrays
        for (std::size_t i : {2u, 3u, 4u, 5u})
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("MATLAB:mxVTKWriter:invalidInputs",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a double array.");
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

#endif // MEX_FUNCTION

//----------------------------------------------------------------------------------------
