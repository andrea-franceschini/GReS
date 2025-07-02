#include "VTUWriter.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

#include "mex.h"
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs < 7 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidInputs", "At least 6 input arguments are required." );
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
  int const nPoints = mxGetM( mxCoord );
  int const dim = mxGetN( mxCoord );

  mxArray const * const mxCells = prhs[3];
  double const * cells = (double *) mxGetPr( mxCells );
  mxArray const * const mxCellVTKType = prhs[4];
  double const * cellVTKType = (double *) mxGetPr( mxCellVTKType );
  mxArray const * const mxCellNumVerts = prhs[5];
  double const * cellNumVerts = (double *) mxGetPr( mxCellNumVerts );

  int const nCells = mxGetM( mxCells );
  int const maxCellNumVerts = mxGetN( mxCells );

  int nConnections = 0;
  for( int i = 0; i < nCells*maxCellNumVerts; ++i )
  {
    nConnections += ( (int)cells[i] != 0 ) ? 1 : 0;
  }

  int * const cellNumVertsInt = (int *) malloc( nCells*sizeof( int ) );
  int * const cellVTKTypeInt = (int *) malloc( nCells*sizeof( int ) );
  for( int i = 0; i < nCells; ++i )
  {
    cellNumVertsInt[i] = (int) cellNumVerts[i];
    cellVTKTypeInt[i] = (int) cellVTKType[i];
  }

  int * const cellsInt = (int *) malloc( nConnections*sizeof( int ) );
  int k = 0;
  for( int i = 0; i < nCells; ++i )
  {
    for( int j = 0; j < cellNumVertsInt[i]; ++j )
    {
      cellsInt[k++] = (int) cells[i+nCells*j] - 1;
    }
  }

  for( int i = 0; i < nCells; ++i )
  {
    if( cellVTKTypeInt[i] != VTK_TRIANGLE &&
       cellVTKTypeInt[i] != VTK_QUAD &&
       cellVTKTypeInt[i] != VTK_TETRA &&
       cellVTKTypeInt[i] != VTK_HEXAHEDRON &&
       cellVTKTypeInt[i] != VTK_HEXAHEDRON2 &&
       cellVTKTypeInt[i] != VTK_QUAD2)
    {
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
          mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:struct", "Scalar and vector data has to be saved as struct." );
        }
        mwSize const nFields = mxGetNumberOfFields( ptr );
        if( nFields != 2 )
        {
          mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:struct", "Struct for scalar and vector data must have 2 fields." );
        }
        mwSize const nStructElems = mxGetNumberOfElements( ptr );
        for( mwIndex i = 0; i < nStructElems; i++ )
        {
          mxArray const * const fieldName = mxGetField( ptr, i, "name" );
          if( fieldName == NULL )
          {
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "name field must exist." );
          }
          if( !mxIsChar( fieldName ) )
          {
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "name field must by a string." );
          }

          char fieldNameC[bufSize];
          status = mxGetString( fieldName, fieldNameC, bufSize );
          if( status != 0 )
          {
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:fieldName", "Not enough space. String is truncated." );
          }

          mxArray const * const data = mxGetField( ptr, i, "data" );
          if( data == NULL )
          {
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "data field must exist." );
          }
          double const * const ptrData = (double *) mxGetData( data );
          int const size = mxGetM( data );
          if( size != nPoints && size != nCells )
          {
            mexErrMsgIdAndTxt( "MATLAB:mxVTKWriter:invalidField", "Size inconsistency." );
          }
          int const dim = mxGetN( data );
          if( size == nPoints )
          {
            vtuFile.add_field< double >( std::string( fieldNameC ), ptrData, size, dim );
          }
          else
          {
            vtuFile.add_cell_field< double >( std::string( fieldNameC ), ptrData, size, dim );
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
