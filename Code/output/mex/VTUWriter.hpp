#ifndef VTU_WRITER_HPP
#define VTU_WRITER_HPP

#include <cstdint>
#include <fstream>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

constexpr int VTK_TETRA       = 10;
constexpr int VTK_TRIANGLE    = 5;
constexpr int VTK_QUAD        = 9;
constexpr int VTK_HEXAHEDRON  = 12;
constexpr int VTK_HEXAHEDRON2 = 29;
constexpr int VTK_QUAD2       = 28;
constexpr int VTK_POLYGON     = 7;
constexpr int VALUES_IN_COLUMN = 10;

enum class VTUFormat
{
  ASCII,
  BINARY
};

template<typename T>
inline char const* vtkTypeName();

template<>
inline char const* vtkTypeName<int>() { return "Int32"; }

template<>
inline char const* vtkTypeName<double>() { return "Float64"; }

class VTUWriter
{
public:
  void set_mesh(int dim,
                int numPoints,
                double const* coord,
                int numCells,
                int numConnections,
                int const* cells,
                int const* cellVTKType,
                int const* nVertices);

  void set_time(double time);
  void set_format(VTUFormat format);

  bool write_mesh(std::string const& path);

  template<typename T>
  void add_field(std::string const& name,
                 T const* data,
                 int size,
                 int dimension)
  {
    addDataNode(m_pointData<T>(), name, data, size, dimension);
  }

  template<typename T>
  void add_cell_field(std::string const& name,
                      T const* data,
                      int size,
                      int dimension)
  {
    addDataNode(m_cellData<T>(), name, data, size, dimension);
  }

  void clear();

private:
  template<typename T>
  class VTKDataNode
  {
  public:
    VTKDataNode(std::string name,
                T const* data,
                int size,
                int nComponents)
      : m_name(std::move(name)),
        m_data(data),
        m_size(size),
        m_numComponents(nComponents)
    {}

    std::uint64_t payloadBytes() const
    {
      return static_cast<std::uint64_t>(m_size) *
             static_cast<std::uint64_t>(m_numComponents) * sizeof(T);
    }

    void writeAsciiHeader(std::ostream& os) const
    {
      os << "<DataArray type=\"" << vtkTypeName<T>()
         << "\" Name=\"" << m_name
         << "\" NumberOfComponents=\"" << m_numComponents
         << "\" format=\"ascii\">\n";
    }

    void writeBinaryHeader(std::ostream& os, std::uint64_t offset) const
    {
      os << "<DataArray type=\"" << vtkTypeName<T>()
         << "\" Name=\"" << m_name
         << "\" NumberOfComponents=\"" << m_numComponents
         << "\" format=\"appended\" offset=\"" << offset << "\"/>\n";
    }

    void writeAsciiPayload(std::ostream& os) const
    {
      if (m_numComponents == 1)
      {
        int column = 0;
        for (int i = 0; i < m_size; ++i)
        {
          os << m_data[i];
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
      }
      else
      {
        for (int i = 0; i < m_size; ++i)
        {
          for (int c = 0; c < m_numComponents; ++c)
          {
            if (c != 0) os << ' ';
            os << m_data[i + c * m_size];
          }
          os << '\n';
        }
      }
    }

    void writeBinaryBlock(std::ostream& os) const
    {
      std::uint64_t const nBytes = payloadBytes();
      os.write(reinterpret_cast<char const*>(&nBytes), sizeof(nBytes));

      if (m_numComponents == 1)
      {
        os.write(reinterpret_cast<char const*>(m_data),
                 static_cast<std::streamsize>(nBytes));
        return;
      }

      std::vector<T> interleaved(
        static_cast<std::size_t>(m_size) * m_numComponents);
      for (int i = 0; i < m_size; ++i)
      {
        for (int c = 0; c < m_numComponents; ++c)
        {
          interleaved[static_cast<std::size_t>(i) * m_numComponents + c] =
            m_data[i + c * m_size];
        }
      }
      os.write(reinterpret_cast<char const*>(interleaved.data()),
               static_cast<std::streamsize>(nBytes));
    }

  private:
    std::string m_name;
    T const* m_data;
    int m_size;
    int m_numComponents;
  };

  template<typename T>
  static void addDataNode(std::vector<VTKDataNode<T>>& nodes,
                          std::string const& name,
                          T const* data,
                          int size,
                          int dimension)
  {
    nodes.emplace_back(name, data, size, dimension);
  }

  template<typename T>
  std::vector<VTKDataNode<T>>& m_pointData();

  template<typename T>
  std::vector<VTKDataNode<T>>& m_cellData();

  void writeAscii(std::ostream& os);
  void writeBinary(std::ostream& os);
  void writeHeader(std::ostream& os) const;
  void writeTimeAscii(std::ostream& os) const;
  void writePointsAscii(std::ostream& os) const;
  void writeCellsAscii(std::ostream& os) const;

  template<typename T>
  static void writeAsciiNodes(std::ostream& os,
                              std::vector<VTKDataNode<T>> const& nodes)
  {
    for (auto const& node : nodes)
    {
      node.writeAsciiHeader(os);
      node.writeAsciiPayload(os);
      os << "</DataArray>\n";
    }
  }

  int m_dim = 0;
  int m_numPoints = 0;
  double const* m_points = nullptr;
  int m_numCells = 0;
  int m_numConnections = 0;
  int const* m_cells = nullptr;
  int const* m_cellVTKType = nullptr;
  int const* m_numVertices = nullptr;
  double m_time = 0.0;
  VTUFormat m_format = VTUFormat::ASCII;

  std::vector<VTKDataNode<int>> m_pointDataInt;
  std::vector<VTKDataNode<int>> m_cellDataInt;
  std::vector<VTKDataNode<double>> m_pointDataDouble;
  std::vector<VTKDataNode<double>> m_cellDataDouble;
};

template<>
inline std::vector<VTUWriter::VTKDataNode<int>>& VTUWriter::m_pointData<int>()
{
  return m_pointDataInt;
}

template<>
inline std::vector<VTUWriter::VTKDataNode<double>>& VTUWriter::m_pointData<double>()
{
  return m_pointDataDouble;
}

template<>
inline std::vector<VTUWriter::VTKDataNode<int>>& VTUWriter::m_cellData<int>()
{
  return m_cellDataInt;
}

template<>
inline std::vector<VTUWriter::VTKDataNode<double>>& VTUWriter::m_cellData<double>()
{
  return m_cellDataDouble;
}

#endif