#include "VTUWriter.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <vector>

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
     << "<DataArray type=\"Float64\" Name=\"TIME\" "
        "NumberOfTuples=\"1\" format=\"ascii\">\n"
     << m_time << '\n'
     << "</DataArray>\n"
     << "</FieldData>\n";
}

void VTUWriter::writePointsAscii(std::ostream& os) const
{
  os << "<Points>\n"
     << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        "format=\"ascii\">\n";

  for (int i = 0; i < m_numPoints; ++i)
  {
    double const x = m_points[i];
    double const y = m_points[i + m_numPoints];
    double const z =
      m_dim == 3 ? m_points[i + 2 * m_numPoints] : 0.0;

    os << x << ' ' << y << ' ' << z << '\n';
  }

  os << "</DataArray>\n"
     << "</Points>\n";
}

void VTUWriter::writeCellsAscii(std::ostream& os) const
{
  os << "<Cells>\n"
     << "<DataArray type=\"Int32\" Name=\"connectivity\" "
        "format=\"ascii\">\n";

  int column = 0;
  for (int i = 0; i < m_numConnections; ++i)
  {
    os << m_cells[i];
    if (++column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else
    {
      os << ' ';
    }
  }
  if (column != 0) os << '\n';
  os << "</DataArray>\n";

  int minType = 0;
  int maxType = 0;
  if (m_numCells > 0)
  {
    auto const minmax =
      std::minmax_element(m_cellVTKType, m_cellVTKType + m_numCells);
    minType = *minmax.first;
    maxType = *minmax.second;
  }

  os << "<DataArray type=\"UInt8\" Name=\"types\" "
        "format=\"ascii\" RangeMin=\""
     << minType << "\" RangeMax=\"" << maxType << "\">\n";

  column = 0;
  for (int i = 0; i < m_numCells; ++i)
  {
    os << m_cellVTKType[i];
    if (++column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else
    {
      os << ' ';
    }
  }
  if (column != 0) os << '\n';
  os << "</DataArray>\n";

  os << "<DataArray type=\"Int64\" Name=\"offsets\" "
        "format=\"ascii\">\n";

  std::int64_t cumulative = 0;
  column = 0;
  for (int i = 0; i < m_numCells; ++i)
  {
    cumulative += static_cast<std::int64_t>(m_numVertices[i]);
    os << cumulative;
    if (++column == VALUES_IN_COLUMN)
    {
      os << '\n';
      column = 0;
    }
    else
    {
      os << ' ';
    }
  }
  if (column != 0) os << '\n';

  os << "</DataArray>\n"
     << "</Cells>\n";
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
    for (auto const& node : m_pointDataInt) node.writeAscii(os);
    for (auto const& node : m_pointDataDouble) node.writeAscii(os);
    os << "</PointData>\n";
  }

  writeCellsAscii(os);

  if (!m_cellDataInt.empty() || !m_cellDataDouble.empty())
  {
    os << "<CellData>\n";
    for (auto const& node : m_cellDataInt) node.writeAscii(os);
    for (auto const& node : m_cellDataDouble) node.writeAscii(os);
    os << "</CellData>\n";
  }

  os << "</Piece>\n"
     << "</UnstructuredGrid>\n"
     << "</VTKFile>\n";
}

// ---------------------------------------------------------------------------
// True binary output: every <DataArray> only carries a format="appended"
// offset="..." reference; the actual bytes live in a single
// <AppendedData encoding="raw"> block at the end of the file, each array
// prefixed by an 8-byte (UInt64) size header. No base64, no text encoding
// of numeric data at all.
// ---------------------------------------------------------------------------
void VTUWriter::writeAppended(std::ostream& os)
{
  // Build the a few derived buffers exactly once, up front, so the
  // descriptor pass and the payload pass can both reference them by
  // pointer without re-deriving anything.
  std::vector<double> points(static_cast<std::size_t>(m_numPoints) * 3);
  for (int i = 0; i < m_numPoints; ++i)
  {
    std::size_t const base = 3 * static_cast<std::size_t>(i);
    points[base]     = m_points[i];
    points[base + 1] = m_points[i + m_numPoints];
    points[base + 2] =
      m_dim == 3 ? m_points[i + 2 * m_numPoints] : 0.0;
  }

  std::vector<std::int32_t> connectivity(
    static_cast<std::size_t>(m_numConnections));
  for (int i = 0; i < m_numConnections; ++i)
  {
    connectivity[static_cast<std::size_t>(i)] =
      static_cast<std::int32_t>(m_cells[i]);
  }

  std::vector<std::uint8_t> types(static_cast<std::size_t>(m_numCells));
  std::vector<std::int64_t> offsets(static_cast<std::size_t>(m_numCells));
  std::int64_t cumulative = 0;
  for (int i = 0; i < m_numCells; ++i)
  {
    types[static_cast<std::size_t>(i)] =
      static_cast<std::uint8_t>(m_cellVTKType[i]);
    cumulative += static_cast<std::int64_t>(m_numVertices[i]);
    offsets[static_cast<std::size_t>(i)] = cumulative;
  }

  // Running offset (relative to the first byte after the '_' marker) and
  // the list of payload writers, collected in the exact order their
  // descriptors are emitted below.
  std::uint64_t currentOffset = 0;
  std::vector<std::function<void(std::ostream&)>> payloadWriters;

  auto emitDataArray = [&os, &currentOffset, &payloadWriters](
                          char const* type,
                          char const* name,
                          int numComponents,
                          char const* extraAttrs,
                          std::uint64_t byteSize,
                          std::function<void(std::ostream&)> writePayload)
  {
    os << "<DataArray type=\"" << type << "\"";
    if (name != nullptr && name[0] != '\0') os << " Name=\"" << name << "\"";
    os << " NumberOfComponents=\"" << numComponents << "\"";
    if (extraAttrs != nullptr && extraAttrs[0] != '\0') os << ' ' << extraAttrs;
    os << " format=\"appended\" offset=\"" << currentOffset << "\"/>\n";

    payloadWriters.push_back(std::move(writePayload));
    currentOffset += sizeof(std::uint64_t) + byteSize;
  };

  // Generic helper: write an 8-byte size header followed by the raw bytes
  // of a contiguous buffer.
  auto writeRaw = [](std::ostream& os2, void const* data, std::uint64_t byteSize)
  {
    os2.write(reinterpret_cast<char const*>(&byteSize), sizeof(byteSize));
    if (byteSize > 0)
      os2.write(reinterpret_cast<char const*>(data),
                static_cast<std::streamsize>(byteSize));
  };

  writeHeader(os);

  // ---- FieldData (TIME) ----
  os << "<FieldData>\n";
  emitDataArray("Float64", "TIME", 1, "NumberOfTuples=\"1\"", sizeof(double),
    [this, writeRaw](std::ostream& os2) { writeRaw(os2, &m_time, sizeof(double)); });
  os << "</FieldData>\n";

  os << "<Piece NumberOfPoints=\"" << m_numPoints
     << "\" NumberOfCells=\"" << m_numCells << "\">\n";

  // ---- Points ----
  os << "<Points>\n";
  {
    std::uint64_t const byteSize =
      static_cast<std::uint64_t>(points.size()) * sizeof(double);
    emitDataArray("Float64", nullptr, 3, "", byteSize,
      [&points, byteSize, writeRaw](std::ostream& os2)
      { writeRaw(os2, points.data(), byteSize); });
  }
  os << "</Points>\n";

  // ---- PointData ----
  if (!m_pointDataInt.empty() || !m_pointDataDouble.empty())
  {
    os << "<PointData>\n";
    for (auto const& node : m_pointDataInt)
    {
      emitDataArray(vtkTypeName<int>(), node.name().c_str(),
        node.numComponents(), "", node.appendedPayloadSize(),
        [&node](std::ostream& os2) { node.writeAppendedPayload(os2); });
    }
    for (auto const& node : m_pointDataDouble)
    {
      emitDataArray(vtkTypeName<double>(), node.name().c_str(),
        node.numComponents(), "", node.appendedPayloadSize(),
        [&node](std::ostream& os2) { node.writeAppendedPayload(os2); });
    }
    os << "</PointData>\n";
  }

  // ---- Cells ----
  os << "<Cells>\n";
  {
    std::uint64_t const byteSize =
      static_cast<std::uint64_t>(connectivity.size()) * sizeof(std::int32_t);
    emitDataArray("Int32", "connectivity", 1, "", byteSize,
      [&connectivity, byteSize, writeRaw](std::ostream& os2)
      { writeRaw(os2, connectivity.data(), byteSize); });
  }
  {
    std::uint64_t const byteSize =
      static_cast<std::uint64_t>(types.size()) * sizeof(std::uint8_t);
    emitDataArray("UInt8", "types", 1, "", byteSize,
      [&types, byteSize, writeRaw](std::ostream& os2)
      { writeRaw(os2, types.data(), byteSize); });
  }
  {
    std::uint64_t const byteSize =
      static_cast<std::uint64_t>(offsets.size()) * sizeof(std::int64_t);
    emitDataArray("Int64", "offsets", 1, "", byteSize,
      [&offsets, byteSize, writeRaw](std::ostream& os2)
      { writeRaw(os2, offsets.data(), byteSize); });
  }
  os << "</Cells>\n";

  // ---- CellData ----
  if (!m_cellDataInt.empty() || !m_cellDataDouble.empty())
  {
    os << "<CellData>\n";
    for (auto const& node : m_cellDataInt)
    {
      emitDataArray(vtkTypeName<int>(), node.name().c_str(),
        node.numComponents(), "", node.appendedPayloadSize(),
        [&node](std::ostream& os2) { node.writeAppendedPayload(os2); });
    }
    for (auto const& node : m_cellDataDouble)
    {
      emitDataArray(vtkTypeName<double>(), node.name().c_str(),
        node.numComponents(), "", node.appendedPayloadSize(),
        [&node](std::ostream& os2) { node.writeAppendedPayload(os2); });
    }
    os << "</CellData>\n";
  }

  os << "</Piece>\n"
     << "</UnstructuredGrid>\n";

  // ---- Raw payloads, in the same order the descriptors were emitted ----
  os << "<AppendedData encoding=\"raw\">\n_";
  for (auto const& writePayload : payloadWriters)
  {
    writePayload(os);
  }
  os << "\n</AppendedData>\n"
     << "</VTKFile>\n";
}

bool VTUWriter::write_mesh(std::string const& path)
{
  std::ofstream os;
  std::vector<char> buffer(4 * 1024 * 1024);
  os.rdbuf()->pubsetbuf(
    buffer.data(), static_cast<std::streamsize>(buffer.size()));

  os.open(path, std::ios::out | std::ios::binary | std::ios::trunc);
  if (!os) return false;

  if (m_format == VTUFormat::BINARY)
  {
    writeAppended(os);
  }
  else
  {
    writeAscii(os);
  }

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