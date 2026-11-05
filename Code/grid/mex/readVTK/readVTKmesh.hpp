#ifndef READ_VTK_MESH
#define READ_VTK_MESH

#include <vector>
#include <string>

// Structure to store point coordinates
struct Point {
    float x, y, z;
};

// Structure to store cell information
struct Cell {
    std::vector<int> pointIndices;
    int cellType;
    int cellTag;
};

// Function to read VTK Unstructured Grid (ASCII)
void readVTKmesh(const std::string& filename, std::vector<Point>& points, std::vector<Cell>& cells);

#endif // READ_VTK_MESH
