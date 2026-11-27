% mode 1 (classic master/slave approach)



g1 = structuredMesh(1,1,1,[0 1],[0 1],[0 1]);
g2 = structuredMesh(2,2,2,[1 2],[0 1],[0 1]);
gI = structuredMesh(0,3,3,1,[0 1],[0 1]);

g1.processGeometry;
g2.processGeometry;

mat = Materials();
mat.addSolid('name',"rock",'cellTags',1);
mat.addConstitutiveLaw("rock","Elastic",'youngModulus',5e3,'poissonRatio',0.0);

bc1 = Boundaries(g1);

bc1.addBC( ...
    'name',           "fix", ...
    'type',           "dirichlet", ...
    'field',          "surface", ...
    'entityListType', "tags", ...
    'entityList',     5, ...
    'components', [1,2,3], ...
    'variable',       "displacements");

bc1.addBCEvent("fix", 'time', 0.0, 'value', 0.0);

bc2 = Boundaries(g2);

bc2.addBC( ...
    'name',           "load", ...
    'type',           "neumann", ...
    'field',          "surface", ...
    'entityListType', "tags", ...
    'entityList',     6, ...
    'components', 1, ...
    'variable',       "displacements");

bc2.addBCEvent("load", 'time', 0.0, 'value', -1.0);


d1 = Domain('grid',g1,'materials',mat,'boundaries',bc1);
d2 = Domain('grid',g2,'materials',mat,'boundaries',bc2);
d1.addPhysicsSolver("Poromechanics");
d2.addPhysicsSolver("Poromechanics");




interf = Interface(gI,[d1,d2]);


interf.grid.findGridIntersections;

interf.addConstraint("MortarContinuity",'quadrature',struct('type',"SegmentBasedQuadrature"))

names = interf.getConstraintNames;
constr = interf.getConstraint("MortarContinuity");


model = Model([d1;d2],interf);

% Use a constant time step over the interval [0, 10].
input = struct( ...
    'Start',           0.0, ...
    'End',             1.0, ...
    'DtInit',          1.0, ...
    'DtMax',           1.0, ...
    'DtMin',           1.0, ...
    'incrementFactor', 1.0);

simParam = SimulationParameters(input);

% Define the requested output times and the output file prefix.

out = OutState( ...
    'printTimes',  1, ...
    'outputFile', "test");

solver = NonLinearImplicit( ...
    'simulationparameters', simParam, ...
    'model',              model, ...
    'output',               out);

solver.simulationLoop();


