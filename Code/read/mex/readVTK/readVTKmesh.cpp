#include "readVTKmesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

void readVTKmesh(const std::string& filename, std::vector<Point>& points, std::vector<Cell>& cells) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Could not open file.");
    }

    std::string line;
    std::getline(file, line); // Read header
    if (line.find("# vtk DataFile Version 2.0") == std::string::npos) {
        throw std::runtime_error("Error: Unsupported VTK file version.");
    }

    std::getline(file, line); // Title (ignored)
    std::getline(file, line); // Format (must be ASCII)
    if (line.find("ASCII") == std::string::npos) {
        throw std::runtime_error("Error: Only ASCII format is supported.");
    }

    std::getline(file, line); // Dataset type
    if (line.find("DATASET UNSTRUCTURED_GRID") == std::string::npos) {
        throw std::runtime_error("Error: Only UNSTRUCTURED_GRID is supported.");
    }

    // Read POINTS
    int numPoints;
    std::string dataType;

    std::getline(file, line); // POINTS 4092 double
    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword >> numPoints >> dataType;
    if (keyword != "POINTS") {
        throw std::runtime_error("Error: Expected POINTS section.");
    }
    
    points.resize(numPoints);
    
    //file >> dataType; // Read data type
    for (int i = 0; i < numPoints; i++) {
        file >> points[i].x >> points[i].y >> points[i].z;
    }

    // Read CELLS
    int numCells, totalIndices;
    file >> line >> numCells >> totalIndices;
    if (line != "CELLS") {
        throw std::runtime_error("Error: Expected CELLS section.");
    }

    cells.resize(numCells);
    for (int i = 0; i < numCells; i++) {
        int numIndices;
        file >> numIndices;
        cells[i].pointIndices.resize(numIndices);
        for (int j = 0; j < numIndices; j++) {
            file >> cells[i].pointIndices[j];
        }
    }

    // Read CELL_TYPES
    int cellTypeCount;
    file >> line >> cellTypeCount;
    if (line != "CELL_TYPES" || cellTypeCount != numCells) {
        throw std::runtime_error("Error: CELL_TYPES section mismatch.");
    }

    for (int i = 0; i < numCells; i++) {
        file >> cells[i].cellType;
    }

    // Read CELL_DATA (if present)
    while (file >> line) {
        if (line == "CELL_DATA") {
            int numCellData;
            file >> numCellData;
            if (numCellData != numCells) {
                throw std::runtime_error("Error: CELL_DATA count mismatch.");
            }

            file >> line; // Read "SCALARS"
            if (line != "SCALARS") {
                throw std::runtime_error("Error: Expected SCALARS keyword in CELL_DATA.");
            }

            std::string dataName, dataType;
            int numComponents;
            file >> dataName >> dataType >> numComponents;
            if (numComponents != 1) {
                throw std::runtime_error("Error: Only 1-component scalars are supported.");
            }

            file >> line; // Read "LOOKUP_TABLE"
            std::string lookupTableName;
            file >> lookupTableName; // Lookup table name

            // Read scalar values
            for (int i = 0; i < numCells; i++) {
                file >> cells[i].cellTag;
            }
        }
    }
}
