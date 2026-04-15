#include "VTUWriter.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

#include "mex.hpp"
#include "mexAdapter.hpp"

#include <memory>
#include <string>
#include <vector>

class MexFunction : public matlab::mex::Function
{
public:
  void operator()( matlab::mex::ArgumentList outputs,
                   matlab::mex::ArgumentList inputs ) override
  {
    if( inputs.size() < 7 )
    {
      error( "MATLAB:mxVTKWriter:invalidInputs", "At least 7 input arguments are required." );
    }
    else if( outputs.size() > 0 )
    {
      error( "MATLAB:mxVTKWriter:invalidOutputs", "Too many output arguments." );
    }

    matlab::data::Array const fileNameArr = inputs[0];
    if( fileNameArr.getType() != matlab::data::ArrayType::CHAR )
    {
      error( "MATLAB:mxVTKWriter:inputNotString", "Input must be a string." );
    }

    std::string const fileName = toAsciiString( fileNameArr );

    matlab::data::TypedArray<double> const timeArr = getDoubleArray( inputs[1], "time must be a real double scalar." );
    if( timeArr.getNumberOfElements() != 1 )
    {
      error( "MATLAB:mxVTKWriter:time", "time must be a scalar." );
    }
    double const time = *timeArr.begin();

    matlab::data::TypedArray<double> const mxCoord = getDoubleArray( inputs[2], "coord must be a real double matrix." );
    auto const coordDims = mxCoord.getDimensions();
    if( coordDims.size() != 2 )
    {
      error( "MATLAB:mxVTKWriter:coord", "coord must be a 2D matrix." );
    }

    double const * coord = getPointer( mxCoord );
    int const nPoints = static_cast<int>( coordDims[0] );
    int const dim = static_cast<int>( coordDims[1] );

    matlab::data::TypedArray<double> const mxCells = getDoubleArray( inputs[3], "cells must be a real double array." );
    double const * cells = getPointer( mxCells );

    matlab::data::TypedArray<double> const mxCellVTKType = getDoubleArray( inputs[4], "cellVTKType must be a real double array." );
    double const * cellVTKType = getPointer( mxCellVTKType );

    matlab::data::TypedArray<double> const mxCellNumVerts = getDoubleArray( inputs[5], "cellNumVerts must be a real double array." );
    double const * cellNumVerts = getPointer( mxCellNumVerts );

    int const nCells = static_cast<int>( mxCellNumVerts.getNumberOfElements() );
    int const nCellTypes = static_cast<int>( mxCellVTKType.getNumberOfElements() );
    int const nConnections = static_cast<int>( mxCells.getNumberOfElements() );

    if( nCellTypes != nCells )
    {
      error( "MATLAB:mxVTKWriter:invalidField", "cellVTKType and cellNumVerts must have the same number of elements." );
    }

    std::vector<int> cellNumVertsInt( nCells );
    std::vector<int> cellVTKTypeInt( nCells );

    int totalConnections = 0;
    for( int i = 0; i < nCells; ++i )
    {
      cellNumVertsInt[i] = static_cast<int>( cellNumVerts[i] );
      cellVTKTypeInt[i] = static_cast<int>( cellVTKType[i] );

      if( cellNumVertsInt[i] <= 0 )
      {
        error( "MATLAB:mxVTKWriter:cellNumVerts", "Each entry of cellNumVerts must be positive." );
      }

      totalConnections += cellNumVertsInt[i];
    }

    if( totalConnections != nConnections )
    {
      error( "MATLAB:mxVTKWriter:connectivity", "Length of flat connectivity must match sum(cellNumVerts)." );
    }

    std::vector<int> cellsInt( nConnections );
    for( int i = 0; i < nConnections; ++i )
    {
      int const id = static_cast<int>( cells[i] );
      if( id <= 0 || id > nPoints )
      {
        error( "MATLAB:mxVTKWriter:connectivity", "Connectivity contains an invalid node index." );
      }
      cellsInt[i] = id - 1;
    }

    for( int i = 0; i < nCells; ++i )
    {
      if( cellVTKTypeInt[i] != VTK_TRIANGLE &&
          cellVTKTypeInt[i] != VTK_QUAD &&
          cellVTKTypeInt[i] != VTK_TETRA &&
          cellVTKTypeInt[i] != VTK_HEXAHEDRON &&
          cellVTKTypeInt[i] != VTK_HEXAHEDRON2 &&
          cellVTKTypeInt[i] != VTK_QUAD2 &&
          cellVTKTypeInt[i] != VTK_POLYGON )
      {
        error( "MATLAB:mxVTKWriter:cellType", "Element type not yet supported." );
      }
    }

    VTUWriter vtuFile;
    vtuFile.set_mesh( dim,
                      nPoints,
                      coord,
                      nCells,
                      nConnections,
                      cellsInt.data(),
                      cellVTKTypeInt.data(),
                      cellNumVertsInt.data() );
    vtuFile.set_time( time );

    for( std::size_t iStruct = 0; iStruct < inputs.size() - 6; ++iStruct )
    {
      matlab::data::Array const ptr = inputs[6 + iStruct];
      if( ptr.getNumberOfElements() > 0 )
      {
        if( ptr.getType() != matlab::data::ArrayType::STRUCT )
        {
          error( "MATLAB:mxVTKWriter:struct", "Scalar and vector data has to be saved as struct." );
        }

        matlab::data::StructArray const s = ptr;
        for( auto const& elem : s )
        {
          matlab::data::Array const fieldName = elem["name"];
          if( fieldName.isEmpty() )
          {
            error( "MATLAB:mxVTKWriter:invalidField", "name field must exist." );
          }

          matlab::data::Array const dataArr = elem["data"];
          if( dataArr.isEmpty() )
          {
            error( "MATLAB:mxVTKWriter:invalidField", "data field must exist." );
          }

          if( fieldName.getType() != matlab::data::ArrayType::CHAR )
          {
            error( "MATLAB:mxVTKWriter:invalidField", "name field must by a string." );
          }

          std::string const fieldNameC = toAsciiString( fieldName );

          matlab::data::TypedArray<double> const data =
              getDoubleArray( dataArr, "data field must be a real double array." );

          auto const dataDims = data.getDimensions();
          if( dataDims.size() != 2 )
          {
            error( "MATLAB:mxVTKWriter:invalidField", "data field must be a 2D matrix." );
          }

          double const * const ptrData = getPointer( data );
          int const size = static_cast<int>( dataDims[0] );
          int const fieldDim = static_cast<int>( dataDims[1] );

          if( size != nPoints && size != nCells )
          {
            error( "MATLAB:mxVTKWriter:invalidField", "Size inconsistency." );
          }

          if( size == nPoints )
          {
            vtuFile.add_field<double>( fieldNameC, ptrData, size, fieldDim );
          }
          else
          {
            vtuFile.add_cell_field<double>( fieldNameC, ptrData, size, fieldDim );
          }
        }
      }
    }

    if( !vtuFile.write_mesh( fileName ) )
    {
      error( "MATLAB:mxVTKWriter:write", "Failed to write output file." );
    }
  }

private:
  std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
  matlab::data::ArrayFactory factory;

  void error( std::string const& id, std::string const& msg )
  {
    matlabPtr->feval( u"error", 0,
      std::vector<matlab::data::Array>{
        factory.createScalar( id ),
        factory.createScalar( msg )
      } );
  }

  matlab::data::TypedArray<double> getDoubleArray( matlab::data::Array const& arr,
                                                   char const* msg )
  {
    if( arr.getType() != matlab::data::ArrayType::DOUBLE )
    {
      error( "MATLAB:mxVTKWriter:type", msg );
    }
    return matlab::data::TypedArray<double>( arr );
  }

  std::string toAsciiString( matlab::data::Array const& arr )
  {
    matlab::data::CharArray const charArr = arr;
    std::u16string const s16 = charArr.toUTF16();
    return std::string( s16.begin(), s16.end() );
  }

  double const* getPointer( matlab::data::TypedArray<double> const& arr )
  {
    return arr.begin().operator->();
  }
};

#endif // MEX_FUNCTION