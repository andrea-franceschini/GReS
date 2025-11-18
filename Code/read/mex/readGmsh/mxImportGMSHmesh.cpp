#include "readGMSHmesh.hpp"

#define MEX_FUNCTION

#ifdef MEX_FUNCTION

#include "mex.h"
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  if( nrhs != 1 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh:invalidInputs", "Only one input is required." );
  }
  else if( nlhs > 3 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh:undefinedOutputs", "Too many output arguments." );
  }
  else if( mxIsChar( prhs[0] ) != 1 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh:inputNotString", "Input must be a string." );
  }

  int status;
  int const bufSize = 1024;
  char fileName[bufSize];

  status = mxGetString( prhs[0], fileName, bufSize );
  if( status != 0 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh:fileName", "Not enough space. String is truncated." );
  }

  std::vector< Region > regions;
  std::vector< Point > coord;
  std::vector< Element > elems;
  status = readGMSHmesh( std::string( fileName ), coord, elems, regions );
  if( status == 1 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh", "Input file does not exist." );
  }
  else if( status == 2 )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh", "Wrong file format." );
  }

  int const nNodes = coord.size();
  mxArray * mxCoord = mxCreateNumericMatrix( (mwSize) nNodes,
                                             (mwSize) 3,
                                             mxDOUBLE_CLASS,
                                             mxREAL );
  if( mxCoord == NULL )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh:createArray", "Not enough space." );
  }
  double * const ptrCoord = (double *) mxGetData( mxCoord );
  for( int i = 0; i < nNodes; ++i )
  {
    ptrCoord[i] = coord[i].x;
    ptrCoord[nNodes+i] = coord[i].y;
    ptrCoord[2*nNodes+i] = coord[i].z;
  }
  coord.resize( 0 );
  plhs[0] = mxCoord;

  int const nElems = elems.size();
  mxArray * mxElems = mxCreateNumericMatrix( (mwSize) nElems,
                                             (mwSize) MAX_NUM_VERTICES+3,
                                             mxINT32_CLASS,
                                             mxREAL );
  if( mxElems == NULL )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh:createArray", "Not enough space." );
  }
  int * const ptrElems = (int *) mxGetData( mxElems );
  for( int i = 0; i < nElems; ++i )
  {
    int j = 0;
    ptrElems[j+i] = elems[i].ID;
    j += nElems;
    ptrElems[j+i] = elems[i].tag;
    j += nElems;
    ptrElems[j+i] = elems[i].n;
    for( int k = 0; k < elems[i].n; ++k )
    {
      j += nElems;
      ptrElems[j+i] = elems[i].v[k];
    }
    for( int k = elems[i].n; k < MAX_NUM_VERTICES; ++k )
    {
      j += nElems;
      ptrElems[j+i] = 0;
    }
  }
  elems.resize( 0 );
  plhs[1] = mxElems;

  int const nFields = regions.size();
  const char * fnames[] = { "dim", "ID", "name" };
  mxArray * mxRegions = mxCreateStructMatrix( (mwSize) nFields,
                                              (mwSize) 1,
                                              (mwSize) 3,
                                              fnames );
  if( mxRegions == NULL )
  {
    mexErrMsgIdAndTxt( "MATLAB:mxImportGMSHmesh:createArray", "Not enough space." );
  }
  for( int i = 0; i < regions.size(); ++i )
  {
    mxArray * val0 = mxCreateNumericMatrix( (mwSize) 1,
                                            (mwSize) 1,
                                            mxINT32_CLASS,
                                            mxREAL );
    int * const ptr0 = (int *) mxGetData( val0 );
    ptr0[0] = regions[i].dim;
    mxSetField( mxRegions, i, "dim", val0 );
    mxArray * val1 = mxCreateNumericMatrix( (mwSize) 1,
                                            (mwSize) 1,
                                            mxINT32_CLASS,
                                            mxREAL );
    int * const ptr1 = (int *) mxGetData( val1 );
    ptr1[0] = regions[i].ID;
    mxSetField( mxRegions, i, "ID", val1 );
    mxSetField( mxRegions, i, "name", mxCreateString( regions[i].name.c_str() ) );
  }
  regions.resize( 0 );
  plhs[2] = mxRegions;

}

#else // MEX_FUNCTION

#include <iostream>
int main( int argc, char **argv )
{
  if( argc != 2 )
  {
    std::cout << "Only one input is required." << std::endl;
    return 1;
  }
  std::string const fileName( argv[1] );
  std::vector< Point > coord;
  std::vector< Element > elems;
  std::vector< Region > regions;
  readGMSHmesh( fileName, coord, elems, regions );

  for( int i = 0; i < coord.size(); ++i )
  {
    std::cout << coord[i].x << " " << coord[i].y << " " << coord[i].z << std::endl;
  }
  for( int i = 0; i < elems.size(); ++i )
  {
    std::cout << elems[i].ID << " " << elems[i].tag << " ";
    for( int j = 0; j < elems[i].n; ++j )
    {
      std::cout << elems[i].v[j] << " ";
    }
    std::cout << std::endl;
  }
  for( int i = 0; i < regions.size(); ++i )
  {
    std::cout << regions[i].dim << " " << regions[i].ID << " " << regions[i].name << " " << std::endl;
  }
  std::cout << "END" << std::endl;
  return 0;
}

#endif // MEX_FUNCTION
