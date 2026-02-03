
close all;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

%% STEP 1) SINGLE PHYSICS - SINGLE PHASE FLOW

fileName = "01_singlePhysics.xml";

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


%% STEP 2) COUPLING WITH MECHANICS


fileName = "02_moreMaterials.xml";

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

%% STEP 3) USING SUBDOMAINS

%% STEP 4) ADDING A NON-CONFORMING INTERFACE

%% STEP 5) CONTACT MECHANICS AND EMBEDDED FRACTURES