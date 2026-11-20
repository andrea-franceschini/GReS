% mode 1 (classic master/slave approach)

g1 = structuredMesh(1,1,1,[0 1],[0 1],[0 1]);
g2 = structuredMesh(2,2,2,[1 2],[0 1],[0 1]);
gI = structuredMesh(3,3,1,[0 1],[0 1]);

g1.processGeometry;
g2.processGeometry;


d1 = Domain('grid',g1);
d2 = Domain('grid',g2);

interf = Interface([d1,d2]);


interf.grid.findGridIntersections;