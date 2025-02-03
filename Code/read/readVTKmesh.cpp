#include <fstream>
#include <sstream>
#include <regex>
#include "readGMSHmesh.hpp"

int readGMSHmesh( std::string const & fileName,
                  std::vector< Point > & coord,
                  std::vector< Element > & elems,
                  std::vector< Region > & regions )
{
  std::ifstream inFile( fileName );
  if ( !inFile.good() )
  {
    return 1;
  }

  std::string line;
  // Version info
  std::getline( inFile, line );
  if ( line != "$MeshFormat" )
  {
    return 2;
  }
  std::getline( inFile, line );
  std::string version;
  std::stringstream( line ) >> version;
  if ( version != "2.2" )
  {
    return 2;
  }

  while ( inFile >> line )
  {
    if ( line == "$PhysicalNames" )
    {
      int nNames;
      inFile >> nNames;
      // empty string
      std::getline( inFile, line );
      regions.resize( nNames );
      for ( int i = 0; i < nNames; ++i )
      {
        std::string name;
        std::getline( inFile, line );
        std::stringstream( line ) >> regions[i].dim >> regions[i].ID >> name;
        name.erase( std::remove( name.begin(),  name.end(), '"'), name.end() );
        name.erase( std::remove( name.begin(),  name.end(), ';'), name.end() );
        regions[i].name = name;
      }
    }
    else if ( line == "$Nodes" )
    {
      int nNodes;
      inFile >> nNodes;
      // empty string
      std::getline( inFile, line );
      coord.resize( nNodes );
      int j;
      for ( int i = 0; i < nNodes; ++i )
      {
        inFile >> j >> coord[i].x >> coord[i].y >> coord[i].z;
      }
    }
    else if ( line == "$Elements" )
    {
      int nElements;
      inFile >> nElements;
      // empty string
      std::getline( inFile, line );

      elems.resize( nElements );

      int j, k, elemNodes;
      int numUnusedElems = 0;
      int ind = 0;
      for ( int i = 0; i < nElements; ++i )
      {
        std::getline( inFile, line );
        std::stringstream str( line );
        int elemIDread;
        str >> j >> elemIDread;
        bool usedElem = true;
        switch ( elemIDread )
        {
          case ( 1 ):
          {
            // nLin++;
            elemNodes = 2;
            break;
          }
          case ( 2 ):
          {
            // nTri++;
            elemNodes = 3;
            break;
          }
          case ( 3 ):
          {
            // nQua++;
            elemNodes = 4;
            break;
          }
          case ( 4 ):
          {
            // nTet++;
            elemNodes = 4;
            break;
          }
          case ( 5 ):
          {
            // nHex++;
            elemNodes = 8;
            break;
          }
          case ( 6 ):
          {
            // nWed++;
            elemNodes = 6;
            break;
          }
          case ( 7 ):
          {
            // nPyr++;
            elemNodes = 5;
            break;
          }
          default:
          {
            usedElem = false;
            numUnusedElems++;
            break;
          }
        }
        if ( usedElem )
        {
          elems[ind].ID = elemIDread;
          elems[ind].n = elemNodes;
          str >> k >> elems[ind].tag >> k;
          for ( j = 0; j < elemNodes; ++j)
          {
            str >> k;
            elems[ind].v[j] = k;
          }
          ind++;
        }
      }
    }
  }
  inFile.close();

  return 0;
}
