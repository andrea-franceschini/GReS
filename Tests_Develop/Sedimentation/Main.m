close all;
% clear;
input_dir = 'Inputs/';
file_SimP = fullfile(input_dir,'simparam.xml');
file_Mat = fullfile(input_dir,'materials.xml');
file_Mesh = fullfile(input_dir,'mesh.xml');
file_Bcs = fullfile(input_dir,'boundaries.xml');
file_Outp = fullfile(input_dir,'output.xml');
file_Solver = fullfile(input_dir,'solver.xml');

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(file_SimP);

% Create an object of the Materials class and read the materials file
mat = Materials(file_Mat);

% Create the Mesh object
topology = Mesh();

% Import mesh data into the Mesh object
% topology.importMesh(file_Mesh);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Create an object of the "Faces" class and process the face properties
faces = Faces(topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Creating boundaries conditions.
bound = Boundaries(file_Bcs,grid);

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility for the solution
printUtils = OutState(topology,file_Outp);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound,...
                     'OutState',printUtils);

domain.addPhysicsSolver(file_Solver);

% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = 1.e5;
% domain.state.data.potential(:) = domain.state.data.pressure+ mat.getFluid().getFluidSpecWeight()*topology.cellCentroid(:,3);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme.
Solver = FCSolver(simParam,domain);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()