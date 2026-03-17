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
// ALL FIXES APPLIED
// -----------------------------------------------------------------------
// [FIX-A] mexErrMsgIdAndTxt() not declared in the pure C++ MEX API.
//         Replaced with throwError() routing through
//         getEngine()->feval(u"error", ..., createCharArray()).
//
// [FIX-C] ArgumentList::size() / operator[] are NOT const-qualified.
//         validateArguments() takes non-const ArgumentList& references.
//
// [FIX-E] factory.createScalar<T>() requires std::is_arithmetic<T>.
//         All string arguments use factory.createCharArray() instead.
//
// [NEW-1] mxGetString() + fixed 512-byte char[] buffer replaced with the
//         C++ MEX API CharArray::toAscii() pattern. This correctly handles
//         any path length and removes the silent truncation risk.
//         MATLAB_STRING inputs are also accepted (both "string" and 'char'
//         array syntax work from the MATLAB caller side).
//
// [NEW-2] mxCreateDoubleMatrix() + mxGetPr() replaced with
//         factory.createArray<double>({rows, cols}).
//         TypedArray<T> does not expose a raw pointer — elements are
//         accessed via flat column-major index [i + j*nRows], which is
//         identical to the pointer arithmetic in the original:
//           coordPtr[i + j*numPoints]  →  coordTable[i + j*numPoints]
//         This makes the indexing intent explicit and avoids any proxy issues.
//
// [NEW-3] std::runtime_error from readVTKmesh() is caught and rethrown as
//         a MATLAB exception via throwError(), ensuring C++ destructors
//         unwind correctly — unlike mexErrMsgIdAndTxt which longjmps past them.
//
// [NEW-4] Input validation checks both argument count AND that the input
//         is a string type (ArrayType::CHAR or ArrayType::MATLAB_STRING).
//
// [NEW-5] cellMatrix zero-initialised via std::fill so padding columns for
//         cells with fewer than maxVerts nodes are guaranteed to be 0,
//         matching the implicit zero-fill of mxCreateDoubleMatrix().
//
// [NEW-6] std::size_t used for all array sizes and loop bounds.
// -----------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "readVTKmesh.hpp"

#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>  // std::fill

//----------------------------------------------------------------------------------------
// Convenience aliases
//----------------------------------------------------------------------------------------
using namespace matlab::data;
using matlab::mex::ArgumentList;

//----------------------------------------------------------------------------------------
// MexFunction
//----------------------------------------------------------------------------------------
class MexFunction : public matlab::mex::Function {

    ArrayFactory factory;

public:

    void operator()(ArgumentList outputs, ArgumentList inputs) override
    {
        // [FIX-C][NEW-4]
        validateArguments(outputs, inputs);

        // -----------------------------------------------------------------------
        // [NEW-1] Extract filename from MATLAB char array.
        //         CharArray::toAscii() returns a std::string — no fixed buffer,
        //         no truncation risk, works with any path length.
        // -----------------------------------------------------------------------
        const CharArray    filenameArr = inputs[0];
        const std::string  filename    = filenameArr.toAscii();

        // -----------------------------------------------------------------------
        // Call the VTK reader — may throw std::runtime_error
        // -----------------------------------------------------------------------
        std::vector<Point> points;
        std::vector<Cell>  cells;

        try {
            readVTKmesh(filename, points, cells);
        } catch (const std::exception& e) {
            // [NEW-3] Re-throw as MATLAB exception — RAII destructors run cleanly
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
        // [NEW-2] Build coordTable: numPoints x 3 double matrix (column-major).
        //
        //   MATLAB stores matrices column-major, so for a numPoints x 3 matrix:
        //     col 0 (x):  indices [0             .. numPoints-1    ]
        //     col 1 (y):  indices [numPoints      .. 2*numPoints-1  ]
        //     col 2 (z):  indices [2*numPoints    .. 3*numPoints-1  ]
        //
        //   Flat index for row i, col j  →  i + j * numPoints
        //   This is identical to the original: coordPtr[i + j*numPoints]
        // -----------------------------------------------------------------------
        TypedArray<double> coordTable =
            factory.createArray<double>({numPoints, 3});

        for (std::size_t i = 0; i < numPoints; ++i) {
            coordTable[i                  ] = points[i].x;  // col 0
            coordTable[i + numPoints      ] = points[i].y;  // col 1
            coordTable[i + 2 * numPoints  ] = points[i].z;  // col 2
        }

        // -----------------------------------------------------------------------
        // [NEW-2][NEW-5] Build cellMatrix: numCells x (3 + maxVerts) matrix.
        //
        //   Column layout (same as original):
        //     col 0          : VTK cell type
        //     col 1          : cell tag
        //     col 2          : number of vertices in this cell
        //     col 3..2+maxV  : vertex indices (1-based); 0 if cell has fewer verts
        //
        //   Zero-initialise via std::fill so padding slots are 0,
        //   matching the implicit zero-fill of mxCreateDoubleMatrix(). [NEW-5]
        // -----------------------------------------------------------------------
        const std::size_t numCols = 3 + maxVerts;

        TypedArray<double> cellMatrix =
            factory.createArray<double>({numCells, numCols});
        std::fill(cellMatrix.begin(), cellMatrix.end(), 0.0);  // [NEW-5]

        for (std::size_t i = 0; i < numCells; ++i) {
            // col 0: VTK cell type
            cellMatrix[i                      ] = static_cast<double>(cells[i].cellType);
            // col 1: cell tag
            cellMatrix[i +       numCells     ] = static_cast<double>(cells[i].cellTag);
            // col 2: number of vertices
            cellMatrix[i + 2   * numCells     ] = static_cast<double>(cells[i].pointIndices.size());
            // cols 3..: vertex indices, 0-based → 1-based
            for (std::size_t j = 0; j < cells[i].pointIndices.size(); ++j)
                cellMatrix[i + (3 + j) * numCells] =
                    static_cast<double>(cells[i].pointIndices[j] + 1);
        }

        // -----------------------------------------------------------------------
        // Return outputs to MATLAB
        // -----------------------------------------------------------------------
        outputs[0] = std::move(coordTable);
        outputs[1] = std::move(cellMatrix);
    }

private:

    // [FIX-C] non-const refs — ArgumentList methods are not const-qualified
    // [NEW-4] Validates both count and type of the filename argument
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

        // Accept both 'char array' and "string" scalar from MATLAB caller
        const ArrayType t = inputs[0].getType();
        if (t != ArrayType::CHAR && t != ArrayType::MATLAB_STRING)
            throwError("VTK:InputError",
                       "Input argument 1 (fileName) must be a char array or string.");
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

//----------------------------------------------------------------------------------------
