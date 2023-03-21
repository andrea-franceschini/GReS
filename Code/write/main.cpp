#include "VTUWriter.hpp"

int main()
{
  VTUWriter v;
  int a = 3;
  v.add_field( "test", &a, 1, 1 );
  return 0;
}
