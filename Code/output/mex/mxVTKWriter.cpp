#include "VTUWriter.hpp"

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <algorithm>
#include <cctype>
#include <memory>
#include <string>
#include <vector>

class MexFunction : public matlab::mex::Function
{
public:
void operator()(matlab::mex::ArgumentList outputs,
matlab::mex::ArgumentList inputs) override
{
if (inputs.size() < 8)
error("MATLAB:mxVTKWriter:invalidInputs",
"At least 8 inputs are required: filename, time, coordinates, "
"cells, cellVTKType, cellNumVerts, format, and data structs.");
if (!outputs.empty())
error("MATLAB:mxVTKWriter:invalidOutputs", "No outputs are returned.");

std::string const fileName = getString(inputs[0], "filename must be text.");
double const time = getScalarDouble(inputs[1], "time must be a real double scalar.");
std::string format = getString(inputs[6], "format must be 'ascii' or 'binary'.");
std::transform(format.begin(), format.end(), format.begin(),
[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
if (format != "ascii" && format != "binary")
error("MATLAB:mxVTKWriter:format", "format must be 'ascii' or 'binary'.");

auto const mxCoord = getDoubleArray(inputs[2], "coord must be a real double matrix.");
auto const coordDims = mxCoord.getDimensions();
if (coordDims.size() != 2 || (coordDims[1] != 2 && coordDims[1] != 3))
error("MATLAB:mxVTKWriter:coord", "coord must be an nPoints-by-2 or nPoints-by-3 matrix.");

int const nPoints = static_cast<int>(coordDims[0]);
int const dim = static_cast<int>(coordDims[1]);
double const* coord = getPointer(mxCoord);

auto const mxCells = getDoubleArray(inputs[3], "cells must be a real double array.");
auto const mxCellVTKType = getDoubleArray(inputs[4], "cellVTKType must be a real double array.");
auto const mxCellNumVerts = getDoubleArray(inputs[5], "cellNumVerts must be a real double array.");

double const* cells = getPointer(mxCells);
double const* cellVTKType = getPointer(mxCellVTKType);
double const* cellNumVerts = getPointer(mxCellNumVerts);

int const nCells = static_cast<int>(mxCellNumVerts.getNumberOfElements());
int const nCellTypes = static_cast<int>(mxCellVTKType.getNumberOfElements());
int const nConnections = static_cast<int>(mxCells.getNumberOfElements());
if (nCellTypes != nCells)
error("MATLAB:mxVTKWriter:invalidField",
"cellVTKType and cellNumVerts must have the same length.");

std::vector<int> cellNumVertsInt(nCells);
std::vector<int> cellVTKTypeInt(nCells);
int totalConnections = 0;
for (int i = 0; i < nCells; ++i)
{
cellNumVertsInt[i] = static_cast<int>(cellNumVerts[i]);
cellVTKTypeInt[i] = static_cast<int>(cellVTKType[i]);
if (cellNumVertsInt[i] <= 0)
error("MATLAB:mxVTKWriter:cellNumVerts",
"Each cellNumVerts entry must be positive.");
totalConnections += cellNumVertsInt[i];

int const type = cellVTKTypeInt[i];
if (type != VTK_TRIANGLE && type != VTK_QUAD && type != VTK_TETRA &&
type != VTK_HEXAHEDRON && type != VTK_HEXAHEDRON2 &&
type != VTK_QUAD2 && type != VTK_POLYGON)
error("MATLAB:mxVTKWriter:cellType", "Element type not supported.");
}
if (totalConnections != nConnections)
error("MATLAB:mxVTKWriter:connectivity",
"numel(cells) must equal sum(cellNumVerts).");

std::vector<int> cellsInt(nConnections);
for (int i = 0; i < nConnections; ++i)
{
int const id = static_cast<int>(cells[i]);
if (id <= 0 || id > nPoints)
error("MATLAB:mxVTKWriter:connectivity",
"Connectivity contains an invalid one-based node index.");
cellsInt[i] = id - 1;
}

VTUWriter writer;
writer.set_mesh(dim, nPoints, coord, nCells, nConnections,
cellsInt.data(), cellVTKTypeInt.data(),
cellNumVertsInt.data());
writer.set_time(time);
writer.set_format(format == "binary" ? VTUFormat::BINARY : VTUFormat::ASCII);

for (std::size_t inputId = 7; inputId < inputs.size(); ++inputId)
{
matlab::data::Array const input = inputs[inputId];
if (input.getNumberOfElements() == 0) continue;
if (input.getType() != matlab::data::ArrayType::STRUCT)
error("MATLAB:mxVTKWriter:struct",
"Point and cell data must be supplied as struct arrays.");

matlab::data::StructArray const fields = input;
for (auto const& field : fields)
{
matlab::data::Array const nameArray = field["name"];
matlab::data::Array const dataArray = field["data"];
if (nameArray.isEmpty() || dataArray.isEmpty())
error("MATLAB:mxVTKWriter:invalidField",
"Each field requires nonempty 'name' and 'data' entries.");

std::string const name = getString(nameArray, "Field name must be text.");
auto const data = getDoubleArray(dataArray,
"Field data must be a real double matrix.");
auto const dims = data.getDimensions();
if (dims.size() != 2)
error("MATLAB:mxVTKWriter:invalidField",
"Field data must be a two-dimensional matrix.");

int const size = static_cast<int>(dims[0]);
int const nComponents = static_cast<int>(dims[1]);
if (size != nPoints && size != nCells)
error("MATLAB:mxVTKWriter:invalidField",
"Field row count must equal nPoints or nCells.");

double const* ptr = getPointer(data);
if (size == nPoints)
writer.add_field<double>(name, ptr, size, nComponents);
else
writer.add_cell_field<double>(name, ptr, size, nComponents);
}
}

if (!writer.write_mesh(fileName))
error("MATLAB:mxVTKWriter:write", "Failed to write the VTU file.");
}

private:
std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
matlab::data::ArrayFactory factory;

void error(std::string const& id, std::string const& message)
{
matlabPtr->feval(u"error", 0,
std::vector<matlab::data::Array>{factory.createScalar(id),
factory.createScalar(message)});
}

matlab::data::TypedArray<double>
getDoubleArray(matlab::data::Array const& array, char const* message)
{
if (array.getType() != matlab::data::ArrayType::DOUBLE)
{
error("MATLAB:mxVTKWriter:type", message);
}

return matlab::data::TypedArray<double>(array);
}

double getScalarDouble(matlab::data::Array const& array, char const* message)
{
auto const data = getDoubleArray(array, message);
if (data.getNumberOfElements() != 1)
error("MATLAB:mxVTKWriter:scalar", message);
return *data.begin();
}

std::string getString(matlab::data::Array const& array, char const* message)
{
if (array.getType() == matlab::data::ArrayType::CHAR)
{
matlab::data::CharArray const chars = array;
std::u16string const s16 = chars.toUTF16();
return std::string(s16.begin(), s16.end());
}
error("MATLAB:mxVTKWriter:text", message);
return {};
}

double const* getPointer(matlab::data::TypedArray<double> const& array)
{
return array.begin().operator->();
}
};
