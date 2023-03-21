#include <vector>
#include <string>

#define MAX_NUM_VERTICES 8

struct Point
{
  double x, y, z;
};
struct Element
{
  int ID, tag;
  int n;
  int v[MAX_NUM_VERTICES];
};
struct Region
{
  int dim, ID;
  std::string name;
};

int readGMSHmesh( std::string const & fileName,
                  std::vector< Point > & coord,
                  std::vector< Element > & elems,
                  std::vector< Region > & regions );
