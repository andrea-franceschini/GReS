#include "VTUWriter.hpp"

#include <algorithm>
#include <cassert>
#include <limits>

namespace
{
constexpr std::uint64_t BLOCK_HEADER_BYTES = sizeof(std::uint64_t);

inline void advanceOffset(std::uint64_t& offset, std::uint64_t payloadBytes)
{
  offset += BLOCK_HEADER_BYTES + payloadBytes;
}

template<typename T>
void writeRawBlock(std::ostream& os, T const* data, std::uint64_t count)
{
  std::uint64_t const nBytes = count * sizeof(T);
  os.write(reinterpret_cast<char const*>(&nBytes), sizeof(nBytes));
  if (nBytes != 0)
  {
    os.write(reinterpret_cast<char const*>(data),
             static_cast<std::streamsize>(nBytes));
  }
}
}

void VTUWriter::set_mesh(int const dim,
                         int const numPoints,
                         double const* const coord,
                         int const numCells,
                         int const numConnections,
                         int const* const cells,
                         int const* const cellVTKType,
                         int const* const nVertices)
{
  assert(dim == 2 || dim == 3);
  m_dim = dim;
  m_numPoints = numPoints;
  m_points = coord;
  m_numCells = numCells;
  m_numConnections = numConnections;
  m_cells = cells;
  m_cellVTKType = cellVTKType;
  m_numVertices = nVertices;
}

void VTUWriter::set_time(double const time)
{
  m_time = time;
}

void VTUWriter::set_format(VTUFormat const format)
{
  m_format = format;
}

void VTUWriter::writeHeader(std::ostream& os) const
{
  os << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\""
        " byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
     << "<UnstructuredGrid>\n";
}

void VTUWriter::writeTimeAscii(std::ostream& os) const
{
  os << "<FieldData>\n"
     << "<DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" "
        "format=\"ascii\">\n"
     << m_time << "\n"
     << "</DataArray>\n"
     << "</FieldData>\n";
}

void VTUWriter::writePointsAscii(std::ostream& os) const
{
  os << "<Points>\n"
     << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  for (int i = 0; i < m_numPoints; ++i)
  {
    os << m_points[i] << ' ' << m_points[i + m_numPoints] << ' ';
    os << (m_dim == 3 ? m_points[i + 2 * m_numPoints] : 0.0) << '\n';
  }

  os << "</DataArray>\n</Points>\n";
}

void VTUWriter::writeCellsAscii(std::ostream& os) const
{
  os << "<Cells>\n";

  os << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  int column = 0;
  for (int i = 0; i < m_numConnections; ++i)
  {
    os << m_cells[i];
    if (++column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else os << ' ';
  }
  if (column != 0) os << '\n';
  os << "</DataArray>\n";

  int minType = 0;
  int maxType = 0;
  if (m_numCells > 0)
  {
    auto const minmax = std::minmax_element(m_cellVTKType,
                                             m_cellVTKType + m_numCells);
    minType = *minmax.first;
    maxType = *minmax.second;
  }

  os << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" "
     << "RangeMin=\"" << minType << "\" RangeMax=\"" << maxType << "\">\n";
  column = 0;
  for (int i = 0; i < m_numCells; ++i)
  {
    os << m_cellVTKType[i];
    if (++column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else os << ' ';
  }
  if (column != 0) os << '\n';
  os << "</DataArray>\n";

  os << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  std::int64_t offset = 0;
  column = 0;
  for (int i = 0; i < m_numCells; ++i)
  {
    offset += m_numVertices[i];
    os << offset;
    if (++column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else os << ' ';
  }
  if (column != 0) os << '\n';
  os << "</DataArray>\n</Cells>\n";
}

void VTUWriter::writeAscii(std::ostream& os)
{
  writeHeader(os);
  writeTimeAscii(os);
  os << "<Piece NumberOfPoints=\"" << m_numPoints
     << "\" NumberOfCells=\"" << m_numCells << "\">\n";

  writePointsAscii(os);

  if (!m_pointDataInt.empty() || !m_pointDataDouble.empty())
  {
    os << "<PointData>\n";
    writeAsciiNodes(os, m_pointDataInt);
    writeAsciiNodes(os, m_pointDataDouble);
    os << "</PointData>\n";
  }

  writeCellsAscii(os);

  if (!m_cellDataInt.empty() || !m_cellDataDouble.empty())
  {
    os << "<CellData>\n";
    writeAsciiNodes(os, m_cellDataInt);
    writeAsciiNodes(os, m_cellDataDouble);
    os << "</CellData>\n";
  }

  os << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

void VTUWriter::writeBinary(std::ostream& os)
{
  std::vector<double> points(static_cast<std::size_t>(m_numPoints) * 3);
  for (int i = 0; i < m_numPoints; ++i)
  {
    points[3 * static_cast<std::size_t>(i)] = m_points[i];
    points[3 * static_cast<std::size_t>(i) + 1] = m_points[i + m_numPoints];
    points[3 * static_cast<std::size_t>(i) + 2] =
      (m_dim == 3 ? m_points[i + 2 * m_numPoints] : 0.0);
  }

  std::vector<std::uint8_t> types(static_cast<std::size_t>(m_numCells));
  std::vector<std::int64_t> offsets(static_cast<std::size_t>(m_numCells));
  std::int64_t cumulative = 0;
  for (int i = 0; i < m_numCells; ++i)
  {
    types[i] = static_cast<std::uint8_t>(m_cellVTKType[i]);
    cumulative += m_numVertices[i];
    offsets[i] = cumulative;
  }

  std::uint64_t offset = 0;
  std::uint64_t const timeOffset = offset;
  advanceOffset(offset, sizeof(double));
  std::uint64_t const pointsOffset = offset;
  advanceOffset(offset, points.size() * sizeof(double));

  std::vector<std::uint64_t> pointIntOffsets;
  std::vector<std::uint64_t> pointDoubleOffsets;
  for (auto const& node : m_pointDataInt)
  {
    pointIntOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }
  for (auto const& node : m_pointDataDouble)
  {
    pointDoubleOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }

  std::uint64_t const connectivityOffset = offset;
  advanceOffset(offset,
    static_cast<std::uint64_t>(m_numConnections) * sizeof(int));
  std::uint64_t const typesOffset = offset;
  advanceOffset(offset, types.size() * sizeof(std::uint8_t));
  std::uint64_t const offsetsOffset = offset;
  advanceOffset(offset, offsets.size() * sizeof(std::int64_t));

  std::vector<std::uint64_t> cellIntOffsets;
  std::vector<std::uint64_t> cellDoubleOffsets;
  for (auto const& node : m_cellDataInt)
  {
    cellIntOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }
  for (auto const& node : m_cellDataDouble)
  {
    cellDoubleOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }

  writeHeader(os);
  os << "<FieldData>\n"
     << "<DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" "
        "format=\"appended\" offset=\"" << timeOffset << "\"/>\n"
     << "</FieldData>\n"
     << "<Piece NumberOfPoints=\"" << m_numPoints
     << "\" NumberOfCells=\"" << m_numCells << "\">\n"
     << "<Points>\n"
     << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        "format=\"appended\" offset=\"" << pointsOffset << "\"/>\n"
     << "</Points>\n";

  if (!m_pointDataInt.empty() || !m_pointDataDouble.empty())
  {
    os << "<PointData>\n";
    for (std::size_t i = 0; i < m_pointDataInt.size(); ++i)
      m_pointDataInt[i].writeBinaryHeader(os, pointIntOffsets[i]);
    for (std::size_t i = 0; i < m_pointDataDouble.size(); ++i)
      m_pointDataDouble[i].writeBinaryHeader(os, pointDoubleOffsets[i]);
    os << "</PointData>\n";
  }

  os << "<Cells>\n"
     << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""
     << connectivityOffset << "\"/>\n"
     << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""
     << typesOffset << "\"/>\n"
     << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\""
     << offsetsOffset << "\"/>\n"
     << "</Cells>\n";

  if (!m_cellDataInt.empty() || !m_cellDataDouble.empty())
  {
    os << "<CellData>\n";
    for (std::size_t i = 0; i < m_cellDataInt.size(); ++i)
      m_cellDataInt[i].writeBinaryHeader(os, cellIntOffsets[i]);
    for (std::size_t i = 0; i < m_cellDataDouble.size(); ++i)
      m_cellDataDouble[i].writeBinaryHeader(os, cellDoubleOffsets[i]);
    os << "</CellData>\n";
  }

  os << "</Piece>\n</UnstructuredGrid>\n"
     << "<AppendedData encoding=\"raw\">\n_";

  writeRawBlock(os, &m_time, 1);
  writeRawBlock(os, points.data(), points.size());
  for (auto const& node : m_pointDataInt) node.writeBinaryBlock(os);
  for (auto const& node : m_pointDataDouble) node.writeBinaryBlock(os);
  writeRawBlock(os, m_cells, static_cast<std::uint64_t>(m_numConnections));
  writeRawBlock(os, types.data(), types.size());
  writeRawBlock(os, offsets.data(), offsets.size());
  for (auto const& node : m_cellDataInt) node.writeBinaryBlock(os);
  for (auto const& node : m_cellDataDouble) node.writeBinaryBlock(os);

  os << "\n</AppendedData>\n</VTKFile>\n";
}

bool VTUWriter::write_mesh(std::string const& path)
{
  std::ofstream os;
  std::vector<char> buffer(4 * 1024 * 1024);
  os.rdbuf()->pubsetbuf(buffer.data(),
                        static_cast<std::streamsize>(buffer.size()));
  os.open(path, std::ios::out | std::ios::binary);
  if (!os) return false;

  if (m_format == VTUFormat::BINARY) writeBinary(os);
  else writeAscii(os);

  os.flush();
  bool const success = static_cast<bool>(os);
  clear();
  return success;
}

void VTUWriter::clear()
{
  m_pointDataInt.clear();
  m_cellDataInt.clear();
  m_pointDataDouble.clear();
  m_cellDataDouble.clear();
}