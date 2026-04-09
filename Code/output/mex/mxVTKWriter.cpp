#include "VTUWriter.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

#include "mex.h"
#include <cstdlib>

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs < 7 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidInputs", "At least 7 input arguments are required." );
  }
  else if( nlhs > 0 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidOutputs", "Too many output arguments." );
  }
  else if( mxIsChar( prhs[0] ) != 1 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:inputNotString", "Input must be a string." );
  }

  int status;
  int const bufSize = 1024;
  char fileName[bufSize];

  status = mxGetString( prhs[0], fileName, bufSize );
  if( status != 0 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:fileName", "Not enough space. String is truncated." );
  }

  double const time = mxGetScalar( prhs[1] );

  mxArray const * const mxCoord = prhs[2];
  double * coord = (double *) mxGetPr( mxCoord );
  int const nPoints = (int) mxGetM( mxCoord );
  int const dim = (int) mxGetN( mxCoord );

  // prhs[3] = flat connectivity vector
  // prhs[4] = cell VTK type
  // prhs[5] = number of vertices for each cell

  mxArray const * const mxCells = prhs[3];
  double const * cells = (double *) mxGetPr( mxCells );

  mxArray const * const mxCellVTKType = prhs[4];
  double const * cellVTKType = (double *) mxGetPr( mxCellVTKType );

  mxArray const * const mxCellNumVerts = prhs[5];
  double const * cellNumVerts = (double *) mxGetPr( mxCellNumVerts );

  int const nCells = (int) mxGetNumberOfElements( mxCellNumVerts );
  int const nCellTypes = (int) mxGetNumberOfElements( mxCellVTKType );
  int const nConnections = (int) mxGetNumberOfElements( mxCells );

  if( nCellTypes != nCells )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "cellVTKType and cellNumVerts must have the same number of elements." );
  }

  int * const cellNumVertsInt = (int *) malloc( nCells*sizeof( int ) );
  int * const cellVTKTypeInt = (int *) malloc( nCells*sizeof( int ) );

  if( cellNumVertsInt == NULL || cellVTKTypeInt == NULL )
  {
    free( cellNumVertsInt );
    free( cellVTKTypeInt );
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:allocation", "Memory allocation failed." );
  }

  int totalConnections = 0;
  for( int i = 0; i < nCells; ++i )
  {
    cellNumVertsInt[i] = (int) cellNumVerts[i];
    cellVTKTypeInt[i] = (int) cellVTKType[i];

    if( cellNumVertsInt[i] <= 0 )
    {
      free( cellNumVertsInt );
      free( cellVTKTypeInt );
      mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:cellNumVerts", "Each entry of cellNumVerts must be positive." );
    }

    totalConnections += cellNumVertsInt[i];
  }

  if( totalConnections != nConnections )
  {
    free( cellNumVertsInt );
    free( cellVTKTypeInt );
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:connectivity", "Length of flat connectivity must match sum(cellNumVerts)." );
  }

  int * const cellsInt = (int *) malloc( nConnections*sizeof( int ) );
  if( cellsInt == NULL )
  {
    free( cellNumVertsInt );
    free( cellVTKTypeInt );
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:allocation", "Memory allocation failed." );
  }

  for( int i = 0; i < nConnections; ++i )
  {
    int const id = (int) cells[i];
    if( id <= 0 || id > nPoints )
    {
      free( cellNumVertsInt );
      free( cellVTKTypeInt );
      free( cellsInt );
      mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:connectivity", "Connectivity contains an invalid node index." );
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
        cellVTKTypeInt[i] != VTK_POLYGON
      )
    {
      free( cellNumVertsInt );
      free( cellVTKTypeInt );
      free( cellsInt );
      mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:cellType", "Element type not yet supported." );
    }
  }

  VTUWriter vtuFile;
  vtuFile.set_mesh( dim, nPoints, coord, nCells, nConnections, cellsInt, cellVTKTypeInt, cellNumVertsInt );
  vtuFile.set_time( time );

  if( nrhs > 6 )
  {
    for( int iStruct = 0; iStruct < ( nrhs-6 ); ++iStruct )
    {
      mxArray const * const ptr = prhs[6+iStruct];
      if( mxGetNumberOfElements( ptr ) > 0 )
      {
        if( mxGetClassID( ptr ) != mxSTRUCT_CLASS )
        {
          free( cellNumVertsInt );
          free( cellVTKTypeInt );
          free( cellsInt );
          mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:struct", "Scalar and vector data has to be saved as struct." );
        }
        mwSize const nFields = mxGetNumberOfFields( ptr );
        if( nFields != 2 )
        {
          free( cellNumVertsInt );
          free( cellVTKTypeInt );
          free( cellsInt );
          mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:struct", "Struct for scalar and vector data must have 2 fields." );
        }
        mwSize const nStructElems = mxGetNumberOfElements( ptr );
        for( mwIndex i = 0; i < nStructElems; i++ )
        {
          mxArray const * const fieldName = mxGetField( ptr, i, "name" );
          if( fieldName == NULL )
          {
            free( cellNumVertsInt );
            free( cellVTKTypeInt );
            free( cellsInt );
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "name field must exist." );
          }
          if( !mxIsChar( fieldName ) )
          {
            free( cellNumVertsInt );
            free( cellVTKTypeInt );
            free( cellsInt );
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "name field must by a string." );
          }

          char fieldNameC[bufSize];
          status = mxGetString( fieldName, fieldNameC, bufSize );
          if( status != 0 )
          {
            free( cellNumVertsInt );
            free( cellVTKTypeInt );
            free( cellsInt );
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:fieldName", "Not enough space. String is truncated." );
          }

          mxArray const * const data = mxGetField( ptr, i, "data" );
          if( data == NULL )
          {
            free( cellNumVertsInt );
            free( cellVTKTypeInt );
            free( cellsInt );
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "data field must exist." );
          }
          double const * const ptrData = (double *) mxGetData( data );
          int const size = (int) mxGetM( data );
          if( size != nPoints && size != nCells )
          {
            free( cellNumVertsInt );
            free( cellVTKTypeInt );
            free( cellsInt );
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "Size inconsistency." );
          }
          int const fieldDim = (int) mxGetN( data );
          if( size == nPoints )
          {
            vtuFile.add_field< double >( std::string( fieldNameC ), ptrData, size, fieldDim );
          }
          else
          {
            vtuFile.add_cell_field< double >( std::string( fieldNameC ), ptrData, size, fieldDim );
          }
        }
      }
    }
  }

  vtuFile.write_mesh( fileName );

  free( cellNumVertsInt );
  free( cellVTKTypeInt );
  free( cellsInt );
}

#endif // MEX_FUNCTION
//----------------------------------------------------------------------------------------
