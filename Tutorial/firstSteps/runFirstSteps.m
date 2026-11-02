%% Unit cube with prescribed pressure gradient
% slightly compressible single phase flow model

% unit cube with 4x4x8 hexahedra
grid = structuredMesh(4,4,8,[0 1],[0 1],[0 1]);

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
bc.addBCEvent("topPressure",'time',0.0,'value',10.0);
%bc.addBCEvent("topPressure",'time',10.0,'value',10.0);

% discretizer
domain = Discretizer('Boundaries', bc, ...
                     'Materials',  mat, ...
                     'Grid',       grid);

domain.addPhysicsSolver("SinglePhaseFlowFVTPFA",'targetRegions',1);

input = struct('Start',0.0,'End',100.0,'DtInit',1e-1,'DtMax',1e1,'DtMin',1e-1,'incrementFactor',1.05);
simparams = SimulationParameters(input);

out = OutState('printTimes',0:1:100,'outputFile',"Output/resultsKV",'matFileName',"Output/results");

solver = NonLinearImplicit('simulationparameters',simparams,'domains',domain, 'output', out);
gresLog().setVerbosity(2);
solver.simulationLoop();

%% equivalent model with xml file

clear
clc

params = readInput('firstSteps.xml');
grid = Grid.create(params.Domain.Geometry);
mat = Materials(params.Domain.Materials);
bc = Boundaries(grid, params.Domain.BoundaryConditions); % example for Boundary Conditions
domain = Discretizer('Boundaries', bc, ...
                     'Materials',  mat, ...
                     'Grid',       grid);
domain.addPhysicsSolvers(params.Domain.Solver)
simparams = SimulationParameters(params.SimulationParameters);
output = OutState(params.Output);
solver = NonLinearImplicit('simulationparameters',simparams,'domains',domain, 'output', output);
solver.simulationLoop();

%% using buildModel utility
params = readInput('firstSteps.xml');
domains = buildModel(params);
simparams = SimulationParameters(params.SimulationParameters);
out = OutState(params.Output);
solver = NonLinearImplicit('simulationparameters',simparams,'domains',domain, 'output', out);
solver.simulationLoop();