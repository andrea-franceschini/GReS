//----------------------------------------------------------------------------------------
// readVTKmesh_wrap.cpp
// Modernized MEX gateway — MathWorks C++ MEX API (R2018a+)
// Uses: mex.hpp + mexAdapter.hpp (matlab::data API)
//
// MATLAB signature:
//   [coordTable, cellMatrix] = readVTKmesh(fileName)
//
// Build command:
//   mex -R2018a readVTKmesh_wrap.cpp
//
// FIXES APPLIED (vs previous version)
// -----------------------------------------------------------------------
// [FIX-G] TypedArray flat-index subscripting causes runtime error
//         "Not enough indices provided" on multi-dimensional arrays.
//         When TypedArray is created with {rows, cols}, operator[](i)
//         does NOT perform flat/linear indexing — it expects
//         multi-dimensional subscripts instead.
//
//         All previous direct writes to coordTable[i + j*rows] and
//         cellMatrix[i + j*numCells] are REPLACED with the established
//         pattern used throughout this codebase:
//           1. Fill a flat std::vector<double> with column-major data
//           2. factory.createArray<double>({rows, cols})
//           3. std::copy(vec.begin(), vec.end(), arr.begin())
//         This is correct because TypedArray iterators DO traverse
//         elements in flat column-major order.
//
// All previously documented fixes (FIX-A/B/C/E, NEW-1..5) are retained.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "readVTKmesh.hpp"

#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>  // std::copy, std::fill

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // Extract filename
        // -----------------------------------------------------------------------
        const CharArray   filenameArr = inputs[0];
        const std::string filename    = filenameArr.toAscii();

        // -----------------------------------------------------------------------
        // Call the VTK reader
        // -----------------------------------------------------------------------
        std::vector<Point> points;
        std::vector<Cell>  cells;

        try {
            readVTKmesh(filename, points, cells);
        } catch (const std::exception& e) {
            throwError("VTK:FileError", e.what());
        }

        // -----------------------------------------------------------------------
        // Compute output dimensions
        // -----------------------------------------------------------------------
        const std::size_t numPoints = points.size();
        const std::size_t numCells  = cells.size();

        std::size_t maxVerts = 0;
        for (const auto& cell : cells)
            if (cell.pointIndices.size() > maxVerts)
                maxVerts = cell.pointIndices.size();

        // -----------------------------------------------------------------------
        // [FIX-G] Build coordTable: numPoints x 3
        //
        //   Fill a flat std::vector column-major (col 0 = x, col 1 = y, col 2 = z)
        //   then std::copy into TypedArray — avoids flat-index subscript error.
        //
        //   Column-major layout for {numPoints, 3}:
        //     indices [0             .. numPoints-1    ] → x values (col 0)
        //     indices [numPoints     .. 2*numPoints-1  ] → y values (col 1)
        //     indices [2*numPoints   .. 3*numPoints-1  ] → z values (col 2)
        // -----------------------------------------------------------------------
        std::vector<double> coordVec(numPoints * 3, 0.0);
        for (std::size_t i = 0; i < numPoints; ++i) {
            coordVec[i                  ] = points[i].x;
            coordVec[i +     numPoints  ] = points[i].y;
            coordVec[i + 2 * numPoints  ] = points[i].z;
        }

        TypedArray<double> coordTable = factory.createArray<double>({numPoints, 3});
        std::copy(coordVec.begin(), coordVec.end(), coordTable.begin());

        // -----------------------------------------------------------------------
        // [FIX-G] Build cellMatrix: numCells x (3 + maxVerts)
        //
        //   Column layout (column-major):
        //     col 0          : VTK cell type
        //     col 1          : cell tag
        //     col 2          : number of vertices
        //     col 3..2+maxV  : vertex indices (1-based), 0-padded
        //
        //   std::fill ensures all padding slots are 0 before writing.
        // -----------------------------------------------------------------------
        const std::size_t numCols = 3 + maxVerts;
        std::vector<double> cellVec(numCells * numCols, 0.0);   // zero-filled

        for (std::size_t i = 0; i < numCells; ++i) {
            cellVec[i                    ] = static_cast<double>(cells[i].cellType);
            cellVec[i +       numCells   ] = static_cast<double>(cells[i].cellTag);
            cellVec[i + 2   * numCells   ] = static_cast<double>(cells[i].pointIndices.size());
            for (std::size_t j = 0; j < cells[i].pointIndices.size(); ++j)
                cellVec[i + (3 + j) * numCells] =
                    static_cast<double>(cells[i].pointIndices[j] + 1); // 1-based
        }

        TypedArray<double> cellMatrix =
            factory.createArray<double>({numCells, numCols});
        std::copy(cellVec.begin(), cellVec.end(), cellMatrix.begin());

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(coordTable);
        outputs[1] = std::move(cellMatrix);
    }

private:

    void validateArguments(ArgumentList& outputs, ArgumentList& inputs)
    {
        if (inputs.size() != 1)
            throwError("VTK:InputError",
                       "Expected 1 input argument (fileName), got " +
                       std::to_string(inputs.size()) + ".");

        if (outputs.size() > 2)
            throwError("VTK:OutputError",
                       "Expected at most 2 output arguments, got " +
                       std::to_string(outputs.size()) + ".");

        const ArrayType t = inputs[0].getType();
        if (t != ArrayType::CHAR && t != ArrayType::MATLAB_STRING)
            throwError("VTK:InputError",
                       "Input argument 1 (fileName) must be a char array or string.");
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
