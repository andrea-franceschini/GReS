clear all 
mesh.coordinates = [ -2	-1	-3;
                      2	-1	-3;
                      2	 1	-3;
                     -2	 1	-3;
                     -2	-1	 0;
                      2	-1   0;
                      2   1	 0;
                     -2   1	 0;
                     -2	-1	 3;
                      2	-1	 3;
                      2	 1	 3;
                     -2	 1	 3];
 mesh.cells = [1   2   3   4   5   6   7   8;
               5   6   7   8   9  10  11  12];
 mesh.nCells = 2;
 mesh.nNodes = 12;
 faceTopology = Faces(mesh);