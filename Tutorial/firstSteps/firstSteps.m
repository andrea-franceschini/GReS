%% Unit cube with prescribed pressure gradient
% slightly compressible single phase flow model

% unit cube with 10x10x10 hexahedra
mesh = structuredMesh(10,10,10,[0 1],[0 1],[0 1]);

numbGaussPoints = 2;
grid = struct('topology', mesh, 'cells', Elements(mesh,numbGaussPoints), 'faces', Faces(mesh));

% materials
mat = Materials();
mat.addFluid('dynamicViscosity',1e-3,'specificWeight',9.81,'compressibility',4.4e-7);
mat.addSolid('name',"sand",'cellTags',1);
mat.addPorousRock("sand","permeability",1e-12,"porosity",0.375);
%mat.addConstitutiveLaw("sand","Elastic",'youngModulus',1e3,'poissonRatio',0.25);
%mat.addCapillaryCurves("sand","type","mualem","beta",2.0,"n",1.0,"kappa",1.0);


% boundary conditions
bc = Boundaries(grid);
bc.addBC('name',"bottomPressure",...
        'type',"dirichlet",...
        'field',"surface",...
        'entityListType',"tags",...
        'entityList',1,...
        'variable',"pressure");


% top surfaces have tag 2 in structuredMesh
bc.addBC('name',"topPressure",...
        'type',"dirichlet",...
        'field',"surface",...
        'entityListType',"tags",...
        'entityList',2,...
        'variable',"pressure");

% bottom pressure - single event -> constant 
bc.addBCEvent("bottomPressure",'time',0.0,'value',0.0);

% top pressure - linear variation from t = 0.0 to t = 1.0
bc.addBCEvent("topPressure",'time',0.0,'value',0.0);
bc.addBCEvent("topPressure",'time',10.0,'value',10.0);

% discretizer
domain = Discretizer('Boundaries', bc, ...
                     'Materials',  mat, ...
                     'Grid',       grid);


input = struct('Start',0.0,'End',10.0,'DtInit',1e-1,'DtMax',1e0,'DtMin',1e-1,'incrementFactor',1.1);
simparams = SimulationParameters(input);

out = OutState('printTimes',[1.0,5.0,10.0],'outputFile',"Output/results",'matFileName',"Output/results");

solver = NonLinearImplicit('simulationparameters',simparams,'domain',domain, 'output', out);
gresLog().setVerbosity(2);
solver.simulationLoop();

%% equivalent model with xml file

clear
clc

params = readInput('firstSteps.xml');
mesh = Mesh.create(params.Geometry);
grid = struct('topology', mesh, 'cells', Elements(mesh,nGP), 'faces', Faces(mesh));
mat = Materials(params.Materials);
bc = Boundaries(grid, params.Domain(1).BoundaryConditions); % example for Boundary Conditions
domain = Discretizer('Boundaries', bound, ...
                     'OutState',   printUtils, ...
                     'Materials',  mat, ...
                     'Grid',       grid);
domain.addPhysicsSolvers(params.Solver)
simparams = SimulationParameters(params.SimulationParameters);
output = OutState(params.Output);
solver = NonLinearImplicit('simulationparameters',simparams,'domain',domain, 'output', printUtils);
solver.simulationLoop();