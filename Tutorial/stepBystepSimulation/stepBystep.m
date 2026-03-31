
close all;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

gresLog().setVerbosity(2)

%% STEP 1) SINGLE PHYSICS - SINGLE PHASE FLOW
clear
clc

fileName = "model01/01_singlePhysics.xml";
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);


mesh = structuredMesh(20,20,10,[0 100],[0 100],[0 10]);

gaussOrder = 2;
elems = Elements(mesh,gaussOrder);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);


bound = Boundaries(grid,params.BoundaryConditions);

printUtils = OutState(params.Output);

domain = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain.addPhysicsSolvers(params.Solver);

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();


%% STEP 2) ADDING SUBDOMAINS
clear
clc

fileName = "model02/02_moreMaterials.xml";
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);


mesh = structuredMesh(20,20,10,[0 100],[0 100],[0 10]);


gaussOrder = 2;
elems = Elements(mesh,gaussOrder);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);

% we manually modify the tag for some clay overburden and underburden
over = mesh.cellCentroid(:,3) < 3;
under = mesh.cellCentroid(:,3) > 7;
mesh.cellTag(under) = 2;
mesh.cellTag(over) = 3;

bound = Boundaries(grid,params.BoundaryConditions);

printUtils = OutState(params.Output);

domain = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain.addPhysicsSolvers(params.Solver);

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();


%% STEP 3) COUPLING WITH MECHANICS
clear
clc

fileName = "model03/03_coupledPoromechanics.xml";
params = readInput(fileName);

% Set parameters of the simulation
simParam = SimulationParameters(params.SimulationParameters);

% Create an object of the Materials class and read the materials file
mat = Materials(params.Materials);

mesh = structuredMesh(15,15,10,[0 100],[0 100],[0 10]);

gaussOrder = 2;
elems = Elements(mesh,gaussOrder);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);

% we manually modify the tag for some clay overburden and underburden
over = mesh.cellCentroid(:,3) < 3;
under = mesh.cellCentroid(:,3) > 7;
mesh.cellTag(under) = 2;
mesh.cellTag(over) = 3;

bound = Boundaries(grid,params.BoundaryConditions);

printUtils = OutState(params.Output);

domain = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain.addPhysicsSolvers(params.Solver);

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();



%% STEP 4) ADDING A NON-CONFORMING INTERFACE
clc
clear

%%%%%%%%%%% domain 1

fileName = "model04/04_domRight.xml";
params = readInput(fileName);

mesh1 = structuredMesh(15,15,10,[20 100],[0 100],[0 10]);
elems = Elements(mesh1,2);
faces = Faces(mesh1);
grid = struct('topology',mesh1,'cells',elems,'faces',faces);

% we manually modify the tag for some clay overburden and underburden
over = mesh1.cellCentroid(:,3) < 3;
under = mesh1.cellCentroid(:,3) > 7;
mesh1.cellTag(under) = 2;
mesh1.cellTag(over) = 3;

mat = Materials(params.Materials);
bound = Boundaries(grid,params.BoundaryConditions);

domain1 = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain1.addPhysicsSolvers(params.Solver);


%%%%%%%%%%% domain 2

fileName = "model04/04_domLeft.xml";
params = readInput(fileName);

mesh2 = structuredMesh(6,6,4,[0 20],[0 100],[0 10]);
elems = Elements(mesh2,2);
faces = Faces(mesh2);
grid = struct('topology',mesh2,'cells',elems,'faces',faces);

mat = Materials(params.Materials);
bound = Boundaries(grid,params.BoundaryConditions);

domain2 = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain2.addPhysicsSolvers(params.Solver);

domains = [domain1;domain2];


params = readInput("model04/04_nonConformingInterface.xml");
interfaces = InterfaceSolver.addInterfaces(domains,params.Interface);

%%%%

simParam = SimulationParameters(params.SimulationParameters);

printUtils = OutState(params.Output);

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domains,...
                           'output',printUtils,...
                           'interfaces',interfaces);
solver.simulationLoop();



%% STEP 5) ADDING A FAULT
clc
clear

%%%%%%%%%%% domain 1

fileName = "model05/05_domRight.xml";
params = readInput(fileName);

mesh1 = structuredMesh(15,15,10,[20 100],[0 100],[0 10]);
elems = Elements(mesh1,2);
faces = Faces(mesh1);
grid = struct('topology',mesh1,'cells',elems,'faces',faces);

% we manually modify the tag for some clay overburden and underburden
over = mesh1.cellCentroid(:,3) < 3;
under = mesh1.cellCentroid(:,3) > 7;
mesh1.cellTag(under) = 2;
mesh1.cellTag(over) = 3;

mat = Materials(params.Materials);
bound = Boundaries(grid,params.BoundaryConditions);

domain1 = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain1.addPhysicsSolvers(params.Solver);


%%%%%%%%%%% domain 2

fileName = "model05/05_domLeft.xml";
params = readInput(fileName);

mesh2 = structuredMesh(6,6,4,[0 20],[0 100],[0 10]);
elems = Elements(mesh2,2);
faces = Faces(mesh2);
grid = struct('topology',mesh2,'cells',elems,'faces',faces);

mat = Materials(params.Materials);
bound = Boundaries(grid,params.BoundaryConditions);

domain2 = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain2.addPhysicsSolvers(params.Solver);

domains = [domain1;domain2];


params = readInput("model05/05_fault.xml");

interfaces = InterfaceSolver.addInterfaces(domains,params.Interface);

% set initial traction
% apply initial traction to the interface
tIni = -1.0;
interfaces{1}.state.traction(1:3:end) = tIni;
interfaces{1}.state.iniTraction(1:3:end) = tIni;
interfaces{1}.stateOld.iniTraction(1:3:end) = tIni;
interfaces{1}.stateOld.traction(1:3:end) = tIni;

%%%%

simParam = SimulationParameters(params.SimulationParameters);

printUtils = OutState(params.Output);

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domains,...
                           'output',printUtils,...
                           'interfaces',interfaces);
solver.simulationLoop();

