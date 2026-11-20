% mode 1 (classic master/slave approach)



g1 = structuredMesh(1,1,1,[0 1],[0 1],[0 1]);
g2 = structuredMesh(2,2,2,[1 2],[0 1],[0 1]);
gI = structuredMesh(0,3,3,1,[0 1],[0 1]);

g1.processGeometry;
g2.processGeometry;

mat = Materials();
mat.addSolid('name',"rock",'cellTags',1);
mat.addConstitutiveLaw("rock","Elastic",'youngModulus',5e3,'poissonRatio',0.0);


d1 = Domain('grid',g1);
d2 = Domain('grid',g2);

interf = Interface(gI,[d1,d2]);


interf.grid.findGridIntersections;


%% example patch test for new mortar capabilities

