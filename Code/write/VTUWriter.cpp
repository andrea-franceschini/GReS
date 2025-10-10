#include "VTUWriter.hpp"
#include <cassert>
#include <iostream>
#include <cmath>

void VTUWriter::write_point_data(std::ostream &os)
{
  if( !m_scalar_point_data_present && !m_vector_point_data_present )
  {
    return;
  }
  os << "<PointData>\n";
  for( auto it = m_point_data_int.begin(); it != m_point_data_int.end(); ++it )
  {
    it->write(os);
  }
  for( auto it = m_point_data_dbl.begin(); it != m_point_data_dbl.end(); ++it )
  {
    it->write(os);
  }
  os << "</PointData>\n";
}

void VTUWriter::write_header( const int n_vertices,
                              const int n_elements,
                              std::ostream &os )
{
  os << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
  os << "<UnstructuredGrid>\n";
  write_time( os );
  os << "<Piece NumberOfPoints=\"" << n_vertices << "\" NumberOfCells=\""
     << n_elements << "\">\n";
}

void VTUWriter::write_time( std::ostream &os )
{
  os << "<FieldData>\n";
  os << "<DataArray type=\"Float64\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\" RangeMin=\""
     << m_time << "\" RangeMax=\"" << m_time << "\">\n";
  os << m_time << "\n";
  os << "</DataArray>\n";
  os << "</FieldData>\n";
}

void VTUWriter::write_footer( std::ostream &os )
{
  os << "</Piece>\n";
  os << "</UnstructuredGrid>\n";
  os << "</VTKFile>\n";
}

void VTUWriter::write_points( int const dim,
                              int const numPoints,
                              double const * const points,
                              std::ostream &os )
{
  os << "<Points>\n";
  os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        "format=\"ascii\">\n";

  int k = 0;
  for( int i = 0; i < numPoints; ++i )
  {
    for( int j = 0; j < dim; ++j )
    {
      os << points[i+j*numPoints];
      if( ( ( k+1 ) / VALUES_IN_COLUMN ) * VALUES_IN_COLUMN == k+1 )
      {
        os << "\n";
      }
      else
      {
        os << " ";
      }
      k++;
    }
  }

  os << "</DataArray>\n";
  os << "</Points>\n";
}

void VTUWriter::write_cells( int const numCells,
                             int const numConnections,
                             int const * const cells,
                             int const * const cellVTKType,
                             int const * const n_vertices,
                             std::ostream &os )
{
  os << "<Cells>\n";
  /////////////////////////////////////////////////////////////////////////////
  // List vertex id's i=0, ..., n_vertices associated with each cell c
  os << "<DataArray type=\"Int64\" Name=\"connectivity\" "
        "format=\"ascii\">\n";
  for (int c = 0; c < numConnections; ++c)
  {
    os << cells[c];
    if( ( ( c+1 ) / VALUES_IN_COLUMN ) * VALUES_IN_COLUMN == c+1 )
    {
      os << "\n";
    }
    else
    {
      os << " ";
    }
  }

  int maxVal = cellVTKType[0];
  int minVal = cellVTKType[0];
  for( int i = 1; i < numCells; ++i )
  {
    if( cellVTKType[i] > maxVal )
    {
      maxVal = cellVTKType[i];
    }
    if( cellVTKType[i] < minVal )
    {
      minVal = cellVTKType[i];
    }
  }

  os << "</DataArray>\n";
  /////////////////////////////////////////////////////////////////////////////
  // List the VTK cell type for each mesh element.
  // This assumes a uniform cell type the entire mesh; to generalize, pass
  // or compute the number of vertices per cell and recompute the cell type
  os << "<DataArray type=\"Int8\" Name=\"types\" format=\"ascii\" "
        "RangeMin=\"" << minVal << "\" RangeMax=\"" << maxVal << "\">\n";
  for (int i = 0; i < numCells; ++i)
  {
    os << cellVTKType[i];
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

  /////////////////////////////////////////////////////////////////////////////
  // List offsets to access the vertex indices of the ith cell. Non-trivial
  // if the mesh is a general polyognal mesh.
  os << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" "
        "RangeMin=\""
     << n_vertices[0] << "\" RangeMax=\"" << numConnections << "\">\n";

  // fixed bug here
  int acc = 0;
  for (int i = 0; i < numCells; ++i)
  {
    acc += n_vertices[i];
    os << acc;
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
  /////////////////////////////////////////////////////////////////////////////
  os << "</Cells>\n";
}

void VTUWriter::write_cell_data(std::ostream &os)
{
  if( !m_scalar_cell_data_present && !m_vector_cell_data_present )
  {
    return;
  }

  os << "<CellData>\n";
  for( auto it = m_cell_data_int.begin(); it != m_cell_data_int.end(); ++it )
  {
    it->write(os);
  }
  for( auto it = m_cell_data_dbl.begin(); it != m_cell_data_dbl.end(); ++it )
  {
    it->write(os);
  }
  os << "</CellData>\n";
}

void VTUWriter::clear()
{
  m_point_data_int.clear();
  m_cell_data_int.clear();
  m_point_data_dbl.clear();
  m_cell_data_dbl.clear();
}

void VTUWriter::set_time( double const time )
{
  m_time = time;
}

void VTUWriter::set_mesh( int const dim,
                          int const numPoints,
                          double const * const points,
                          int const numCells,
                          int const numConnections,
                          int const * const cells,
                          int const * const cellVTKType,
                          int const * const n_vertices )
{
  assert( dim > 1 );

  m_dim = dim;
  m_numPoints = numPoints;
  m_points = points;
  m_numCells = numCells;
  m_numConnections = numConnections;
  m_cells = cells;
  m_cellVTKType = cellVTKType;
  m_numVertices = n_vertices;
}


bool VTUWriter::write_mesh( const std::string & path )
{
  std::ofstream os;
  os.open(path.c_str());
  if( !os.good() )
  {
    os.close();
    return false;
  }

  write_header( m_numPoints, m_numCells, os );
  write_points( m_dim, m_numPoints, m_points, os );
  write_point_data( os );
  write_cells( m_numCells, m_numConnections, m_cells, m_cellVTKType, m_numVertices, os );
  write_cell_data( os );
  write_footer( os );
  clear();

  os.close();
  return true;
}

bool VTUWriter::write_point_cloud( std::ostream &os,
                                   int const numPoints,
                                   int const dim,
                                   double const * const points )
{
  write_header( numPoints, 0, os );
  write_points( dim, numPoints, points, os );
  write_point_data( os );
  os << "<Cells>\n";
  /////////////////////////////////////////////////////////////////////////////
  // List vertex id's i=0, ..., n_vertices associated with each cell c
  os << "<DataArray type=\"Int64\" Name=\"connectivity\" "
        "format=\"ascii\">\n";

  os << "</DataArray>\n";
  os << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" "
        "RangeMin=\""
     <<-1e+299 << "\" RangeMax=\"" << 1e+299 << "\">\n";
  os << "</DataArray>\n";
  os << "</Cells>\n";
  write_footer( os );
  clear();
  return true;
}

bool VTUWriter::write_point_cloud( std::string const & path,
                                   int const numPoints,
                                   int const dim,
                                   double const * const points )
{
  std::ofstream os;
  os.open(path.c_str());
  if (!os.good()) {
    os.close();
    return false;
  }
  write_point_cloud( os, numPoints, dim, points );
  os.close();
  return true;
}
