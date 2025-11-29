% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

simParam = SimulationParameters("simParam.xml");
% Create the Mesh object
topology = Mesh();

% Import the mesh data into the Mesh object
topology.importGMSHmesh('Mesh/Column_hexa.msh');

% Create an object of the Materials class and read the materials file
mat = Materials('materials.xml');


gaussOrder = 2;
elems = Elements(topology,gaussOrder);
faces = Faces(topology);

grid = struct('topology',topology,'cells',elems,'faces',faces);

printUtils = OutState(topology,'output.xml');

% Create an object of the "Boundaries" class 
bound = Boundaries('boundaryConditions.xml',grid);

% Create object handling construction of Jacobian and rhs of the model

% create the Discretizer (key-value pair input)
domain = Discretizer('grid',grid,...
                     'materials',mat,...
                     'boundaries',bound,...
                     'outstate',printUtils);

domain.addPhysicsSolver(solverFile);

% manually apply initial conditions
state = domain.getState();

applyTerzaghiIC(state,mat,topology,-10);

solv = FCSolver(simParam,domain);

%calling analytical solution script
%Terzaghi_analytical(topology, mat, 10)

[simState] = solv.NonLinearLoop();
