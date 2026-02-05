
close all;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

gresLog().setVerbosity(2)

%% STEP 1) SINGLE PHYSICS - SINGLE PHASE FLOW

fileName = "model01/01_singlePhysics.xml";

% Set parameters of the simulation
simParam = SimulationParameters(fileName);

% Create an object of the Materials class and read the materials file
mat = Materials(fileName);


mesh = structuredMesh(20,20,10,[0 100],[0 100],[0 10]);

gaussOrder = 2;
elems = Elements(mesh,gaussOrder);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);


bound = Boundaries(fileName,grid);

printUtils = OutState(mesh,fileName);

domain = Discretizer('Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);
domain.addPhysicsSolver(fileName);


Solver = GeneralSolver(simParam,domain);

Solver.NonLinearLoop();

domain.outstate.finalize()


%% STEP 2) ADDING SUBDOMAINS


fileName = "model02/02_moreMaterials.xml";

% Set parameters of the simulation
simParam = SimulationParameters(fileName);

% Create an object of the Materials class and read the materials file
mat = Materials(fileName);


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

bound = Boundaries(fileName,grid);

printUtils = OutState(mesh,fileName);

domain = Discretizer('Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);
domain.addPhysicsSolver(fileName);


Solver = GeneralSolver(simParam,domain);

Solver.NonLinearLoop();

domain.outstate.finalize()

%% STEP 3) COUPLING WITH MECHANICS



fileName = "model03/03_coupledPoromechanics.xml";

% Set parameters of the simulation
simParam = SimulationParameters(fileName);

% Create an object of the Materials class and read the materials file
mat = Materials(fileName);

mesh = structuredMesh(25,25,15,[0 100],[0 100],[0 10]);

gaussOrder = 2;
elems = Elements(mesh,gaussOrder);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);

% we manually modify the tag for some clay overburden and underburden
over = mesh.cellCentroid(:,3) < 3;
under = mesh.cellCentroid(:,3) > 7;
mesh.cellTag(under) = 2;
mesh.cellTag(over) = 3;

bound = Boundaries(fileName,grid);

printUtils = OutState(mesh,fileName);

domain = Discretizer('Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);
domain.addPhysicsSolver(fileName);


Solver = GeneralSolver(simParam,domain);

Solver.NonLinearLoop();

domain.outstate.finalize();


%% STEP 4) ADDING A NON-CONFORMING INTERFACE


fileName = "model04/04_nonConformingInterface.xml";

% Set parameters of the simulation
simParam = SimulationParameters(fileName);
[domain,interfaces] = buildModel(fileName);

% we manually modify the tag for some clay overburden and underburden
over = domain(2).grid.topology.cellCentroid(:,3) < 3;
under = domain(2).grid.topology.cellCentroid(:,3) > 7;
domain(2).grid.topology.cellTag(under) = 2;
domain(2).grid.topology.cellTag(over) = 3;


Solver = GeneralSolver(simParam,domain,interfaces);

Solver.NonLinearLoop();

Solver.finalizeOutput();



%% STEP 5) ADDING A FAULT