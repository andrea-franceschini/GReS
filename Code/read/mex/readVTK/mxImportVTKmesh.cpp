#include "mex.h"
#include "readVTKmesh.hpp"  // Ensure the header is included here
#include <vector>
#include <string>
#include <stdexcept> 

/**
 * MATLAB MEX Entry Point
 * Usage: [coordTable, cellMatrix] = readVTKmesh(fileName)
 */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 1 || !mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("VTK:InputError", "Usage: readVTKmesh(fileName)");
    }

    // Get file name from input
    char filename[512];
    mxGetString(prhs[0], filename, sizeof(filename));

    std::vector<Point> points; 
    std::vector<Cell> cells;

    try {
        readVTKmesh(filename, points, cells);
    } catch (const std::runtime_error& e) {
        mexErrMsgIdAndTxt("VTK:FileError", e.what());
    }

    // Create MATLAB outputs
    size_t numPoints = points.size();
    size_t numCells = cells.size();

    // Create coordTable: numPoints x 3 matrix
    plhs[0] = mxCreateDoubleMatrix(numPoints, 3, mxREAL);
    double* coordPtr = mxGetPr(plhs[0]);
    for (size_t i = 0; i < numPoints; i++) {
        coordPtr[i] = points[i].x;
        coordPtr[i + numPoints] = points[i].y;
        coordPtr[i + 2 * numPoints] = points[i].z;
    }

    // Determine max number of vertices per cell
    size_t maxVerts = 0;
    for (const auto& cell : cells) {
        if (cell.pointIndices.size() > maxVerts)
            maxVerts = cell.pointIndices.size();
    }

    // Create cellMatrix: numCells x (3 + maxVerts) matrix
    plhs[1] = mxCreateDoubleMatrix(numCells, 3 + maxVerts, mxREAL);
    double* cellPtr = mxGetPr(plhs[1]);

    for (size_t i = 0; i < numCells; i++) {
        cellPtr[i] = cells[i].cellType; // VTK Type
        cellPtr[i + numCells] = cells[i].cellTag; // Cell tag
        cellPtr[i + 2 * numCells] = cells[i].pointIndices.size(); // Number of vertices

        for (size_t j = 0; j < cells[i].pointIndices.size(); j++) {
            cellPtr[i + (3 + j) * numCells] = cells[i].pointIndices[j] + 1; // Convert to MATLAB 1-based indexing
        }
    }
}
