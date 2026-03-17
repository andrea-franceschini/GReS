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
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-G] TypedArray flat-index subscripting causes runtime error
//         "Not enough indices provided" on 2D arrays.
//
//         Affected inputs and their fix:
//
//         coord  {nPoints, dim}  — flat idx in coordVec already used ✅
//                                  (was already a vector copy — no change)
//
//         cells  {nCells, maxCellNumVerts}
//                The nConnections loop and the cellsInt build loop both
//                used mxCells[idx] with a computed flat index on the 2D
//                TypedArray — both now replaced with cellsBuf[idx] after
//                copying the TypedArray into std::vector<double> cellsBuf.
//
//         cellVTKType  / cellNumVerts  — these are 1D column vectors so
//                mxCellVTKType[i] / mxCellNumVerts[i] work for a single
//                index, BUT copying to std::vector is safer and consistent.
//
//         field data  {size, dim}  — dataArr was iterated via begin()/end()
//                into dataVec — already correct ✅ (no change)
//
// All previously documented fixes (FIX-A/B/C/E, NEW-1..6) are retained.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "VTUWriter.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <vector>
#include <string>
#include <algorithm>  // std::find, std::copy

using namespace matlab::data;
using matlab::mex::ArgumentList;

static bool hasField(const StructArray& sa, const std::string& name)
{
    for (const auto& fn : sa.getFieldNames())
        if (std::string(fn) == name) return true;
    return false;
}

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Extract filename and time scalar
        // -----------------------------------------------------------------------
        const CharArray   fileNameArr = inputs[0];
        const std::string fileName    = fileNameArr.toAscii();
        const double      time        = static_cast<double>(
                                            TypedArray<double>(inputs[1])[0]);

        // -----------------------------------------------------------------------
        // [FIX-G] Copy ALL multi-dimensional input arrays into std::vector
        //         before any index arithmetic — avoids flat-index runtime error.
        // -----------------------------------------------------------------------

        // coord: {nPoints, dim}
        const TypedArray<double> mxCoord       = inputs[2];
        const auto coordDims                   = mxCoord.getDimensions();
        const int nPoints                      = static_cast<int>(coordDims[0]);
        const int dim                          = static_cast<int>(coordDims[1]);
        std::vector<double> coordBuf(mxCoord.begin(), mxCoord.end());

        // cells: {nCells, maxCellNumVerts}  ← 2D, was causing the crash
        const TypedArray<double> mxCells       = inputs[3];
        const auto cellsDims                   = mxCells.getDimensions();
        const int nCells                       = static_cast<int>(cellsDims[0]);
        const int maxCellNumVerts              = static_cast<int>(cellsDims[1]);
        std::vector<double> cellsBuf(mxCells.begin(), mxCells.end()); // [FIX-G]

        // cellVTKType / cellNumVerts: 1D vectors, copy for consistency
        const TypedArray<double> mxCellVTKType  = inputs[4];
        const TypedArray<double> mxCellNumVerts = inputs[5];
        std::vector<double> cellVTKTypeBuf (mxCellVTKType.begin(),  mxCellVTKType.end());
        std::vector<double> cellNumVertsBuf(mxCellNumVerts.begin(), mxCellNumVerts.end());

        // -----------------------------------------------------------------------
        // Count non-zero connections using cellsBuf [FIX-G]
        // -----------------------------------------------------------------------
        int nConnections = 0;
        {
            const std::size_t total =
                static_cast<std::size_t>(nCells) *
                static_cast<std::size_t>(maxCellNumVerts);
            for (std::size_t idx = 0; idx < total; ++idx)
                if (static_cast<int>(cellsBuf[idx]) != 0) ++nConnections;  // [FIX-G]
        }

        // -----------------------------------------------------------------------
        // Build int conversion arrays
        // -----------------------------------------------------------------------
        std::vector<int> cellNumVertsInt(static_cast<std::size_t>(nCells));
        std::vector<int> cellVTKTypeInt (static_cast<std::size_t>(nCells));

        for (int i = 0; i < nCells; ++i) {
            cellNumVertsInt[static_cast<std::size_t>(i)] =
                static_cast<int>(cellNumVertsBuf[static_cast<std::size_t>(i)]);
            cellVTKTypeInt [static_cast<std::size_t>(i)] =
                static_cast<int>(cellVTKTypeBuf [static_cast<std::size_t>(i)]);
        }

        // Build flat connectivity array, MATLAB 1-based → 0-based, using cellsBuf
        std::vector<int> cellsInt(static_cast<std::size_t>(nConnections));
        {
            int k = 0;
            for (int i = 0; i < nCells; ++i) {
                for (int j = 0; j < cellNumVertsInt[static_cast<std::size_t>(i)]; ++j) {
                    const std::size_t idx =
                        static_cast<std::size_t>(i) +
                        static_cast<std::size_t>(nCells) *
                        static_cast<std::size_t>(j);
                    cellsInt[static_cast<std::size_t>(k++)] =
                        static_cast<int>(cellsBuf[idx]) - 1;  // [FIX-G] use cellsBuf
                }
            }
        }

        // -----------------------------------------------------------------------
        // Validate cell types
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
        // Set up VTUWriter
        // -----------------------------------------------------------------------
        VTUWriter vtuFile;
        vtuFile.set_mesh(dim, nPoints, coordBuf.data(),
                         nCells, nConnections,
                         cellsInt.data(),
                         cellVTKTypeInt.data(),
                         cellNumVertsInt.data());
        vtuFile.set_time(time);

        // -----------------------------------------------------------------------
        // Process optional field structs (inputs 6, 7, ...)
        // -----------------------------------------------------------------------
        const std::size_t nFixed = 6;
        for (std::size_t iStruct = 0; iStruct < inputs.size() - nFixed; ++iStruct) {

            Array inputArg = inputs[nFixed + iStruct];
            if (inputArg.getNumberOfElements() == 0) continue;

            if (inputArg.getType() != ArrayType::STRUCT)
                throwError("MATLAB:mxVTKWriter:struct",
                           "Scalar and vector data must be saved as struct.");

            StructArray sa = inputs[nFixed + iStruct];

            if (sa.getNumberOfFields() != 2)
                throwError("MATLAB:mxVTKWriter:struct",
                           "Struct for scalar and vector data must have exactly 2 fields.");

            if (!hasField(sa, "name"))
                throwError("MATLAB:mxVTKWriter:invalidField",
                           "name field must exist.");
            if (!hasField(sa, "data"))
                throwError("MATLAB:mxVTKWriter:invalidField",
                           "data field must exist.");

            const std::size_t nStructElems = sa.getNumberOfElements();

            for (std::size_t i = 0; i < nStructElems; ++i) {

                const Array nameRaw = sa[i]["name"];
                if (nameRaw.getType() != ArrayType::CHAR &&
                    nameRaw.getType() != ArrayType::MATLAB_STRING)
                    throwError("MATLAB:mxVTKWriter:invalidField",
                               "name field must be a string.");

                const CharArray   nameArr   = sa[i]["name"];
                const std::string fieldName = nameArr.toAscii();

                const Array dataRaw = sa[i]["data"];
                if (dataRaw.getType() != ArrayType::DOUBLE)
                    throwError("MATLAB:mxVTKWriter:invalidField",
                               "data field must be a double array.");

                const TypedArray<double> dataArr = sa[i]["data"];
                const auto dataDims = dataArr.getDimensions();
                const int  dataRows = static_cast<int>(dataDims[0]);
                const int  dataCols = static_cast<int>(dataDims[1]);

                if (dataRows != nPoints && dataRows != nCells)
                    throwError("MATLAB:mxVTKWriter:invalidField",
                               "Size inconsistency in data field.");

                // Copy via iterators — already correct, no FIX-G needed here
                // since we pass the raw pointer of the vector, not index TypedArray
                std::vector<double> dataVec(dataArr.begin(), dataArr.end());

                if (dataRows == nPoints)
                    vtuFile.add_field<double>(fieldName,
                                              dataVec.data(), dataRows, dataCols);
                else
                    vtuFile.add_cell_field<double>(fieldName,
                                                   dataVec.data(), dataRows, dataCols);
            }
        }

        vtuFile.write_mesh(fileName.c_str());
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() < 7)
            throwError("MATLAB:mxVTKWriter:invalidInputs",
                       "At least 7 input arguments are required.");

        if (outputs.size() > 0)
            throwError("MATLAB:mxVTKWriter:invalidOutputs",
                       "Too many output arguments.");

        const ArrayType t = inputs[0].getType();
        if (t != ArrayType::CHAR && t != ArrayType::MATLAB_STRING)
            throwError("MATLAB:mxVTKWriter:inputNotString",
                       "First input (fileName) must be a string or char array.");

        if (inputs[1].getType() != ArrayType::DOUBLE ||
            inputs[1].getNumberOfElements() != 1)
            throwError("MATLAB:mxVTKWriter:invalidInputs",
                       "Second input (time) must be a real double scalar.");

        for (std::size_t i : {2u, 3u, 4u, 5u})
            if (inputs[i].getType() != ArrayType::DOUBLE)
                throwError("MATLAB:mxVTKWriter:invalidInputs",
                           "Input argument " + std::to_string(i + 1) +
                           " must be a double array.");
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

#endif // MEX_FUNCTION

//----------------------------------------------------------------------------------------
