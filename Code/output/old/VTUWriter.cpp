#include "VTUWriter.hpp"

#include <algorithm>
#include <cassert>
#include <limits>

namespace
{

constexpr std::uint64_t BLOCK_HEADER_BYTES = sizeof(std::uint64_t);

inline void advanceOffset(std::uint64_t& offset,
                          std::uint64_t const payloadBytes)
{
  offset += BLOCK_HEADER_BYTES + payloadBytes;
}

template<typename T>
void writeRawBlock(std::ostream& os,
                   T const* const data,
                   std::uint64_t const count)
{
  std::uint64_t const nBytes = count * sizeof(T);

  os.write(
    reinterpret_cast<char const*>(&nBytes),
    sizeof(nBytes));

  if (nBytes > 0)
  {
    os.write(
      reinterpret_cast<char const*>(data),
      static_cast<std::streamsize>(nBytes));
  }
}

} // namespace

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
  os << "<?xml version=\"1.0\"?>\n"
     << "<VTKFile type=\"UnstructuredGrid\" "
        "version=\"1.0\" "
        "byte_order=\"LittleEndian\" "
        "header_type=\"UInt64\">\n"
     << "<UnstructuredGrid>\n";
}

void VTUWriter::writeTimeAscii(std::ostream& os) const
{
  os << "<FieldData>\n"
     << "<DataArray "
        "type=\"Float64\" "
        "Name=\"TIME\" "
        "NumberOfTuples=\"1\" "
        "format=\"ascii\">\n"
     << m_time << '\n'
     << "</DataArray>\n"
     << "</FieldData>\n";
}

void VTUWriter::writePointsAscii(std::ostream& os) const
{
  os << "<Points>\n"
     << "<DataArray "
        "type=\"Float64\" "
        "NumberOfComponents=\"3\" "
        "format=\"ascii\">\n";

  for (int i = 0; i < m_numPoints; ++i)
  {
    double const x = m_points[i];
    double const y = m_points[i + m_numPoints];
    double const z =
      m_dim == 3
        ? m_points[i + 2 * m_numPoints]
        : 0.0;

    os << x << ' '
       << y << ' '
       << z << '\n';
  }

  os << "</DataArray>\n"
     << "</Points>\n";
}

void VTUWriter::writeCellsAscii(std::ostream& os) const
{
  os << "<Cells>\n";

  /*
   * Connectivity.
   */
  os << "<DataArray "
        "type=\"Int32\" "
        "Name=\"connectivity\" "
        "format=\"ascii\">\n";

  int column = 0;

  for (int i = 0; i < m_numConnections; ++i)
  {
    os << m_cells[i];

    ++column;

    if (column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else
    {
      os << ' ';
    }
  }

  if (column != 0)
  {
    os << '\n';
  }

  os << "</DataArray>\n";

  /*
   * VTK cell types.
   */
  int minType = 0;
  int maxType = 0;

  if (m_numCells > 0)
  {
    auto const minmax =
      std::minmax_element(
        m_cellVTKType,
        m_cellVTKType + m_numCells);

    minType = *minmax.first;
    maxType = *minmax.second;
  }

  os << "<DataArray "
        "type=\"UInt8\" "
        "Name=\"types\" "
        "format=\"ascii\" "
        "RangeMin=\"" << minType << "\" "
        "RangeMax=\"" << maxType << "\">\n";

  column = 0;

  for (int i = 0; i < m_numCells; ++i)
  {
    os << m_cellVTKType[i];

    ++column;

    if (column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else
    {
      os << ' ';
    }
  }

  if (column != 0)
  {
    os << '\n';
  }

  os << "</DataArray>\n";

  /*
   * Cell connectivity offsets.
   */
  os << "<DataArray "
        "type=\"Int64\" "
        "Name=\"offsets\" "
        "format=\"ascii\">\n";

  std::int64_t cumulativeOffset = 0;
  column = 0;

  for (int i = 0; i < m_numCells; ++i)
  {
    cumulativeOffset += m_numVertices[i];

    os << cumulativeOffset;

    ++column;

    if (column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else
    {
      os << ' ';
    }
  }

  if (column != 0)
  {
    os << '\n';
  }

  os << "</DataArray>\n"
     << "</Cells>\n";
}

void VTUWriter::writeAscii(std::ostream& os)
{
  writeHeader(os);

  writeTimeAscii(os);

  os << "<Piece "
        "NumberOfPoints=\"" << m_numPoints << "\" "
        "NumberOfCells=\"" << m_numCells << "\">\n";

  writePointsAscii(os);

  if (!m_pointDataInt.empty() ||
      !m_pointDataDouble.empty())
  {
    os << "<PointData>\n";

    writeAsciiNodes(os, m_pointDataInt);
    writeAsciiNodes(os, m_pointDataDouble);

    os << "</PointData>\n";
  }

  writeCellsAscii(os);

  if (!m_cellDataInt.empty() ||
      !m_cellDataDouble.empty())
  {
    os << "<CellData>\n";

    writeAsciiNodes(os, m_cellDataInt);
    writeAsciiNodes(os, m_cellDataDouble);

    os << "</CellData>\n";
  }

  os << "</Piece>\n"
     << "</UnstructuredGrid>\n"
     << "</VTKFile>\n";
}

void VTUWriter::writeBinary(std::ostream& os)
{
  /*
   * MATLAB stores coordinates column-wise:
   *
   * [x_1 ... x_n, y_1 ... y_n, z_1 ... z_n]
   *
   * VTK expects interleaved tuples:
   *
   * [x_1, y_1, z_1, x_2, y_2, z_2, ...]
   */
  std::vector<double> points(
    static_cast<std::size_t>(m_numPoints) * 3);

  for (int i = 0; i < m_numPoints; ++i)
  {
    std::size_t const base =
      3 * static_cast<std::size_t>(i);

    points[base] =
      m_points[i];

    points[base + 1] =
      m_points[i + m_numPoints];

    points[base + 2] =
      m_dim == 3
        ? m_points[i + 2 * m_numPoints]
        : 0.0;
  }

  /*
   * Use explicitly sized integer types so that the binary payload
   * exactly matches the VTK XML type declarations.
   */
  std::vector<std::int32_t> connectivity(
    static_cast<std::size_t>(m_numConnections));

  for (int i = 0; i < m_numConnections; ++i)
  {
    connectivity[static_cast<std::size_t>(i)] =
      static_cast<std::int32_t>(m_cells[i]);
  }

  std::vector<std::uint8_t> types(
    static_cast<std::size_t>(m_numCells));

  std::vector<std::int64_t> cellOffsets(
    static_cast<std::size_t>(m_numCells));

  std::int64_t cumulativeOffset = 0;

  for (int i = 0; i < m_numCells; ++i)
  {
    types[static_cast<std::size_t>(i)] =
      static_cast<std::uint8_t>(m_cellVTKType[i]);

    cumulativeOffset +=
      static_cast<std::int64_t>(m_numVertices[i]);

    cellOffsets[static_cast<std::size_t>(i)] =
      cumulativeOffset;
  }

  /*
   * Compute appended-data offsets.
   *
   * Each block consists of:
   *
   * [UInt64 number of payload bytes][payload]
   *
   * VTK offsets are measured from the byte immediately following
   * the underscore in <AppendedData encoding="raw">_.
   */
  std::uint64_t offset = 0;

  std::uint64_t const timeOffset = offset;

  advanceOffset(
    offset,
    sizeof(double));

  std::uint64_t const pointsOffset = offset;

  advanceOffset(
    offset,
    static_cast<std::uint64_t>(points.size()) *
      sizeof(double));

  std::vector<std::uint64_t> pointIntOffsets;
  pointIntOffsets.reserve(m_pointDataInt.size());

  for (auto const& node : m_pointDataInt)
  {
    pointIntOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }

  std::vector<std::uint64_t> pointDoubleOffsets;
  pointDoubleOffsets.reserve(m_pointDataDouble.size());

  for (auto const& node : m_pointDataDouble)
  {
    pointDoubleOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }

  std::uint64_t const connectivityOffset = offset;

  advanceOffset(
    offset,
    static_cast<std::uint64_t>(connectivity.size()) *
      sizeof(std::int32_t));

  std::uint64_t const typesOffset = offset;

  advanceOffset(
    offset,
    static_cast<std::uint64_t>(types.size()) *
      sizeof(std::uint8_t));

  std::uint64_t const cellOffsetsOffset = offset;

  advanceOffset(
    offset,
    static_cast<std::uint64_t>(cellOffsets.size()) *
      sizeof(std::int64_t));

  std::vector<std::uint64_t> cellIntOffsets;
  cellIntOffsets.reserve(m_cellDataInt.size());

  for (auto const& node : m_cellDataInt)
  {
    cellIntOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }

  std::vector<std::uint64_t> cellDoubleOffsets;
  cellDoubleOffsets.reserve(m_cellDataDouble.size());

  for (auto const& node : m_cellDataDouble)
  {
    cellDoubleOffsets.push_back(offset);
    advanceOffset(offset, node.payloadBytes());
  }

  /*
   * Write the XML description.
   */
  writeHeader(os);

  os << "<FieldData>\n"
     << "<DataArray "
        "type=\"Float64\" "
        "Name=\"TIME\" "
        "NumberOfTuples=\"1\" "
        "format=\"appended\" "
        "offset=\"" << timeOffset << "\"/>\n"
     << "</FieldData>\n";

  os << "<Piece "
        "NumberOfPoints=\"" << m_numPoints << "\" "
        "NumberOfCells=\"" << m_numCells << "\">\n";

  os << "<Points>\n"
     << "<DataArray "
        "type=\"Float64\" "
        "NumberOfComponents=\"3\" "
        "format=\"appended\" "
        "offset=\"" << pointsOffset << "\"/>\n"
     << "</Points>\n";

  if (!m_pointDataInt.empty() ||
      !m_pointDataDouble.empty())
  {
    os << "<PointData>\n";

    for (std::size_t i = 0;
         i < m_pointDataInt.size();
         ++i)
    {
      m_pointDataInt[i].writeBinaryHeader(
        os,
        pointIntOffsets[i]);
    }

    for (std::size_t i = 0;
         i < m_pointDataDouble.size();
         ++i)
    {
      m_pointDataDouble[i].writeBinaryHeader(
        os,
        pointDoubleOffsets[i]);
    }

    os << "</PointData>\n";
  }

  os << "<Cells>\n"
     << "<DataArray "
        "type=\"Int32\" "
        "Name=\"connectivity\" "
        "format=\"appended\" "
        "offset=\"" << connectivityOffset << "\"/>\n"
     << "<DataArray "
        "type=\"UInt8\" "
        "Name=\"types\" "
        "format=\"appended\" "
        "offset=\"" << typesOffset << "\"/>\n"
     << "<DataArray "
        "type=\"Int64\" "
        "Name=\"offsets\" "
        "format=\"appended\" "
        "offset=\"" << cellOffsetsOffset << "\"/>\n"
     << "</Cells>\n";

  if (!m_cellDataInt.empty() ||
      !m_cellDataDouble.empty())
  {
    os << "<CellData>\n";

    for (std::size_t i = 0;
         i < m_cellDataInt.size();
         ++i)
    {
      m_cellDataInt[i].writeBinaryHeader(
        os,
        cellIntOffsets[i]);
    }

    for (std::size_t i = 0;
         i < m_cellDataDouble.size();
         ++i)
    {
      m_cellDataDouble[i].writeBinaryHeader(
        os,
        cellDoubleOffsets[i]);
    }

    os << "</CellData>\n";
  }

  os << "</Piece>\n"
     << "</UnstructuredGrid>\n";

  /*
   * The underscore must immediately follow the opening tag.
   * The first appended binary block starts immediately after it.
   */
  os << "<AppendedData encoding=\"raw\">_";

  /*
   * Write appended blocks in exactly the same order used when
   * computing the offsets above.
   */
  writeRawBlock(
    os,
    &m_time,
    1);

  writeRawBlock(
    os,
    points.data(),
    static_cast<std::uint64_t>(points.size()));

  for (auto const& node : m_pointDataInt)
  {
    node.writeBinaryBlock(os);
  }

  for (auto const& node : m_pointDataDouble)
  {
    node.writeBinaryBlock(os);
  }

  writeRawBlock(
    os,
    connectivity.data(),
    static_cast<std::uint64_t>(connectivity.size()));

  writeRawBlock(
    os,
    types.data(),
    static_cast<std::uint64_t>(types.size()));

  writeRawBlock(
    os,
    cellOffsets.data(),
    static_cast<std::uint64_t>(cellOffsets.size()));

  for (auto const& node : m_cellDataInt)
  {
    node.writeBinaryBlock(os);
  }

  for (auto const& node : m_cellDataDouble)
  {
    node.writeBinaryBlock(os);
  }

  os << "\n"
     << "</AppendedData>\n"
     << "</VTKFile>\n";
}

bool VTUWriter::write_mesh(std::string const& path)
{
  std::ofstream os;

  /*
   * Use a larger stream buffer to reduce write-system-call overhead.
   * The buffer must remain alive until the file stream is closed.
   */
  std::vector<char> buffer(4 * 1024 * 1024);

  os.rdbuf()->pubsetbuf(
    buffer.data(),
    static_cast<std::streamsize>(buffer.size()));

  os.open(
    path,
    std::ios::out |
    std::ios::binary |
    std::ios::trunc);

  if (!os)
  {
    return false;
  }

  if (m_format == VTUFormat::BINARY)
  {
    writeBinary(os);
  }
  else
  {
    writeAscii(os);
  }

  os.flush();

  bool const success =
    static_cast<bool>(os);

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