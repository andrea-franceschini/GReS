
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


grid = structuredMesh(200,20,10,[0 100],[0 100],[0 10]);

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


grid = structuredMesh(20,20,10,[0 100],[0 100],[0 10]);


% we manually modify the tag for some clay overburden and underburden
over = grid.cells.center(:,3) < 3;
under = grid.cells.center(:,3) > 7;
grid.cells.tag(under) = 2;
grid.cells.tag(over) = 3;

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

grid = structuredMesh(15,15,10,[0 100],[0 100],[0 10]);

% we manually modify the tag for some clay overburden and underburden
over = grid.cells.center(:,3) < 3;
under = grid.cells.center(:,3) > 7;
grid.cells.tag(under) = 2;
grid.cells.tag(over) = 3;

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

grid = structuredMesh(15,15,10,[20 100],[0 100],[0 10]);

% we manually modify the tag for some clay overburden and underburden
over = grid.cells.center(:,3) < 3;
under = grid.cells.center(:,3) > 7;
grid.cells.tag(under) = 2;
grid.cells.tag(over) = 3;

mat = Materials(params.Materials);
bound = Boundaries(grid,params.BoundaryConditions);

domain1 = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain1.addPhysicsSolvers(params.Solver);


%%%%%%%%%%% domain 2

fileName = "model04/04_domLeft.xml";
params = readInput(fileName);

grid = structuredMesh(6,6,4,[0 20],[0 100],[0 10]);

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

grid = structuredMesh(15,15,10,[20 100],[0 100],[0 10]);

% we manually modify the tag for some clay overburden and underburden
over = grid.cells.center(:,3) < 3;
under = grid.cells.center(:,3) > 7;
grid.cells.tag(under) = 2;
grid.cells.tag(over) = 3;

mat = Materials(params.Materials);
bound = Boundaries(grid,params.BoundaryConditions);

domain1 = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);
domain1.addPhysicsSolvers(params.Solver);


%%%%%%%%%%% domain 2

fileName = "model05/05_domLeft.xml";
params = readInput(fileName);

grid = structuredMesh(6,6,4,[0 20],[0 100],[0 10]);

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

