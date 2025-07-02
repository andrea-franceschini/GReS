#ifndef VTU_WRITER_HPP
#define VTU_WRITER_HPP

#include <string>
#include <cassert>

#include <fstream>
#include <string>
#include <vector>

#define VALUES_IN_COLUMN 10

static const int VTK_TETRA = 10;
static const int VTK_TRIANGLE = 5;
static const int VTK_QUAD = 9;
static const int VTK_HEXAHEDRON = 12;
static const int VTK_HEXAHEDRON2 = 29;
static const int VTK_QUAD2 = 28;

//static const int VTK_POLYGON = 7;

namespace
{
template< typename T >
inline std::string getType()
{
  return "";
}

template<>
inline std::string getType< int >()
{
  return "Int32";
}

template<>
inline std::string getType< double >()
{
  return "Float64";
}
}

class VTUWriter {
public:

  /**
   * Write volume mesh to an output stream
   *
   * ostream &os                    output stream where to write vtk mesh (ending with .vtu)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. If there are
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size},
   *                                  ...
   *                                  c_{m,1}, c_{m,2},..., c_{m, cell_size}]
   *                                if there are m cells
   *                                (i.e. index c*i corresponds to the ith
   *                                vertex in the cth cell in the mesh
   */
  void set_mesh( int const dim,
                 int const numPoints,
                 double const * const coord,
                 int const numCells,
                 int const numConnections,
                 int const * const cells,
                 int const * const cellVTKType,
                 int const * const n_vertices );

  void set_time( double const time );

  bool write_mesh( std::string const & path );

  /**
   * Write point cloud to a file
   *
   * const string& path             filename to store vtk mesh (ending with .vtu)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. If there are
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size},
   *                                  ...
   *                                  c_{m,1}, c_{m,2},..., c_{m, cell_size}]
   *                                if there are m cells
   *                                (i.e. index c*i corresponds to the ith
   *                                vertex in the cth cell in the mesh
   */
  bool write_point_cloud( std::string const & path,
                          int const dim,
                          int const numPoints,
                          double const * const points );

  /**
   * Write point cloud to a file
   *
   * ostream &os                    output stream where to write vtk mesh (ending with .vtp)
   * const int dim                  ambient dimension (2D or 3D)
   * const int cell_size            number of vertices per cell
   *                                (3 for triangles, 4 for quads and tets, 8
   *                                for hexes)
   * const vector<double>& points   list of point locations. If there are
   *                                n points in the mesh, the format  of the
   *                                vector is:
   *                                  [x_1, y_1, x_2, y_2, ..., x_n, y_n]
   *                                for 2D and
   *                                  [x_1, y_1, z_1, ..., x_n, y_n, z_n]
   *                                for 3D.
   * const vector<int >& elements   list of point indices per cell. Format  of the
   *                                vector is:
   *                                  [c_{1,1}, c_{1,2},..., c_{1, cell_size},
   *                                  ...
   *                                  c_{m,1}, c_{m,2},..., c_{m, cell_size}]
   *                                if there are m cells
   *                                (i.e. index c*i corresponds to the ith
   *                                vertex in the cth cell in the mesh
   */
  bool write_point_cloud( std::ostream &os,
                          int const dim,
                          int const numPoints,
                          double const * const points );

  /**
   * Add a general field to the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be dimension
   *                                values for each point in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension},
   *                                  ...
   *                                  f_{n,1}, f_{n,2},..., f_{n, dimension}]
   *                                if there are n points in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  template< typename T >
  void add_field( std::string const & name,
                             T const * const data,
                             int const size,
                             int const dimension )
  {
    if( dimension == 1 )
    {
      add_scalar_field( name, data, size );
    }
    else
    {
      add_vector_field( name, data, size, dimension );
    }
  }

  /**
   * Add a general cell/element field to the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be dimension
   *                                values for each cell in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension},
   *                                  ...
   *                                  f_{m,1}, f_{m,2},..., f_{m, dimension}]
   *                                if there are m cells in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  template< typename T >
  void add_cell_field( std::string const & name,
                                  T const * const data,
                                  int const size,
                                  int const dimension)
  {
    if( dimension == 1 )
    {
      add_cell_scalar_field< T >( name, data, size );
    }
    else
    {
      add_cell_vector_field< T >( name, data, size, dimension );
    }
  }

  /**
   * Add a scalar field to the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be one
   *                                value for each point in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_1, f_2,..., f_n]
   *                                if there are n points in the mesh
   */
  template< typename T,
            typename std::enable_if< std::is_same< T, int >::value, bool >::type = true >
  void add_scalar_field( std::string const & name,
                                    T const * const data,
                                    int const size )
  {
    m_point_data_int.push_back( VTKDataNode< T >() );
    m_point_data_int.back().initialize( name, getType< T >(), data, size );
    m_scalar_point_data_present = true;
  }

  template< typename T,
            typename std::enable_if< std::is_same< T, double >::value, bool >::type = true >
  void add_scalar_field( std::string const & name,
                                    T const * const data,
                                    int const size )
  {
    m_point_data_dbl.push_back( VTKDataNode< T >() );
    m_point_data_dbl.back().initialize( name, getType< T >(), data, size );
    m_scalar_point_data_present = true;
  }

  /**
   * Add a scalar field to cells/elements of the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be one
   *                                value for each cell in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_1, f_2,..., f_m]
   *                                if there are m cells in the mesh
   */
  template< typename T,
            typename std::enable_if< std::is_same< T, int >::value, bool >::type = true >
  void add_cell_scalar_field( std::string const & name,
                                         T const * const data,
                                         int const size )
  {
    m_cell_data_int.push_back( VTKDataNode< T >() );
    m_cell_data_int.back().initialize( name, getType< T >(), data, size );
    m_scalar_cell_data_present = true;
  }

  template< typename T,
            typename std::enable_if< std::is_same< T, double >::value, bool >::type = true >
  void add_cell_scalar_field( std::string const & name,
                                         T const * const data,
                                         int const size )
  {
    m_cell_data_dbl.push_back( VTKDataNode< T >() );
    m_cell_data_dbl.back().initialize( name, getType< T >(), data, size );
    m_scalar_cell_data_present = true;
  }

  /**
   * Add a vector field to the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be dimension
   *                                values for each point in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension},
   *                                  ...
   *                                  f_{n,1}, f_{n,2},..., f_{n, dimension}]
   *                                if there are n points in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  template< typename T,
            typename std::enable_if< std::is_same< T, int >::value, bool >::type = true >
  void add_vector_field( std::string const & name,
                                    T const * const data,
                                    int const size,
                                    int const dimension )
  {
    m_point_data_int.push_back( VTKDataNode< T >() );
    m_point_data_int.back().initialize( name, getType< T >(), data, size, dimension );
    m_vector_point_data_present = true;
  }

  template< typename T,
            typename std::enable_if< std::is_same< T, double >::value, bool >::type = true >
  void add_vector_field( std::string const & name,
                                    T const * const data,
                                    int const size,
                                    int const dimension )
  {
    m_point_data_dbl.push_back( VTKDataNode< T >() );
    m_point_data_dbl.back().initialize( name, getType< T >(), data, size, dimension );
    m_vector_point_data_present = true;
  }

  /**
   * Add a vector field to cells/elements of the mesh
   * const string& name             name of the field to store vtk mesh
   * const vector<double>& data     list of field values. There must be dimension
   *                                values for each cell in the mesh to be written.
   *                                Format of the vector is
   *                                  [f_{1,1}, f_{1,2},..., f_{1, dimension},
   *                                  ...
   *                                  f_{m,1}, f_{m,2},..., f_{m, dimension}]
   *                                if there are m cells in the mesh
   * const int dimension            ambient dimension (2D or 3D)
   */
  template< typename T,
            typename std::enable_if< std::is_same< T, int >::value, bool >::type = true >
  void add_cell_vector_field( std::string const & name,
                                         T const * const data,
                                         int const size,
                                         int const dimension )
  {
    m_cell_data_int.push_back( VTKDataNode< T >() );
    m_cell_data_int.back().initialize( name, getType< T >(), data, size, dimension );
    m_vector_cell_data_present = true;
  }

  template< typename T,
            typename std::enable_if< std::is_same< T, double >::value, bool >::type = true >
  void add_cell_vector_field( std::string const & name,
                                         T const * const data,
                                         int const size,
                                         int const dimension )
  {
    m_cell_data_dbl.push_back( VTKDataNode< T >() );
    m_cell_data_dbl.back().initialize( name, getType< T >(), data, size, dimension );
    m_vector_cell_data_present = true;
  }

  // Remove all fields and initialized data from the writer.
  void clear();

private:

  template <typename T>
  class VTKDataNode {
  public:
    VTKDataNode() {}

    VTKDataNode( std::string const & name,
                 std::string const & numeric_type,
                 T const * const data = NULL,
                 int const size = 0,
                 int const n_components = 1 ) :
      m_name( name ),
      m_numeric_type( numeric_type ),
      m_data( data ),
      m_size( size ),
      m_numComponents( n_components )
    {}

    inline T const * const data()
    {
      return m_data;
    }

    void initialize( std::string const & name,
                     std::string const & numeric_type,
                     T const * const data,
                     int const size,
                     int const n_components = 1 )
    {
      m_name = name;
      m_numeric_type = numeric_type;
      m_data = data;
      m_size = size;
      m_numComponents = n_components;
    }

    void write(std::ostream &os) const
    {
      if( m_numComponents == 1 )
      {
        os << "<DataArray type=\"" << m_numeric_type << "\" Name=\"" << m_name
           << "\" NumberOfComponents=\"" << m_numComponents << "\" format=\"ascii\">\n";

        for( int i = 0; i < m_size; ++i )
        {
          os << m_data[i];
          if( ( ( i+1 ) / VALUES_IN_COLUMN ) * VALUES_IN_COLUMN == i+1 )
          {
            os << "\n";
          }
          else
          {
            os << " ";
          }
        }

        os << "</DataArray>\n";
      }
      else
      {
        os << "<DataArray type=\"" << m_numeric_type << "\" Name=\"" << m_name
           << "\" NumberOfComponents=\"" << m_numComponents << "\" format=\"ascii\">\n";

        for( int i = 0; i < m_size; ++i )
        {
          for( int j = 0; j < m_numComponents; ++j )
          {
            os << m_data[i+j*m_size];
            if (j < m_numComponents - 1)
            {
              os << " ";
            }
          }
          os << "\n";
        }

        os << "</DataArray>\n";

        for( int j = 0; j < m_numComponents; ++j )
        {
          os << "<DataArray type=\"" << m_numeric_type << "\" Name=\"" << m_name << "_" << j
             << "\" NumberOfComponents=\"" << 1 << "\" format=\"ascii\">\n";
          for( int i = 0; i < m_size; ++i )
          {
            os << m_data[i+j*m_size];
            if( ( ( i+j*m_size+1 ) / VALUES_IN_COLUMN ) * VALUES_IN_COLUMN == i+j*m_size+1 )
            {
              os << "\n";
            }
            else
            {
              os << " ";
            }
          }
          os << "</DataArray>\n";
        }
      }
    }

    inline bool empty() const
    {
      return m_size <= 0;
    }

  private:
    std::string m_name;
    std::string m_numeric_type;
    T const * m_data;
    int m_size;
    int m_numComponents;
  };

  std::vector< VTKDataNode< int > > m_point_data_int;
  std::vector< VTKDataNode< int > > m_cell_data_int;
  std::vector< VTKDataNode< double > > m_point_data_dbl;
  std::vector< VTKDataNode< double > > m_cell_data_dbl;
  bool m_scalar_point_data_present = false;
  bool m_vector_point_data_present = false;
  bool m_scalar_cell_data_present = false;
  bool m_vector_cell_data_present = false;

  int m_dim;
  int m_numPoints;
  double const * m_points;
  int m_numCells;
  int m_numConnections;
  int const * m_cells;
  int const * m_cellVTKType;
  int const * m_numVertices;
  double m_time = 0.0;

  void write_point_data( std::ostream &os );

  void write_cell_data( std::ostream &os );

  void write_header( int const numPoints,
                     int const numCells,
                     std::ostream &os );

  void write_time( std::ostream &os );

  void write_footer( std::ostream &os );

  void write_points( int const numPoints,
                     int const dim,
                     double const * const points,
                     std::ostream &os );

  void write_cells( int const numCells,
                    int const numConnections,
                    int const * const cells,
                    int const * const cellVTKType,
                    int const * const n_vertices,
                    std::ostream &os );

};

#endif // VTU_WRITER_HPP
