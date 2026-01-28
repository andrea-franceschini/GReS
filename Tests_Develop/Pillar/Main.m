close all;
% clear;
output_dir = 'Outputs';
input_dir = 'Inputs';
figures_dir = fullfile(output_dir,"Images");


%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,"Materials",'matTable.xml'));

% Create the Mesh object
topology = Mesh();

% Choosing the mesh file
availMesh = [ "Column.msh", "Column1x1x30.msh", "Column4x4x40.msh","Column22.msh",];
fileName = availMesh(4);

% Import mesh data into the Mesh object
topology.importMesh(fullfile(input_dir,'Mesh',fileName));

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Create an object of the "Faces" class and process the face properties
faces = Faces(topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Creating boundaries conditions.
bound = Boundaries(fullfile(input_dir,'boundaries.xml'),grid);

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility
printUtils = OutState(topology,fullfile(input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound,...
                     'OutState',printUtils);

domain.addPhysicsSolver(fullfile(input_dir,'solver.xml'));

% set initial conditions directly modifying the state object
z = elems.mesh.cellCentroid(:,3);
gamma_w = getFluid(mat).getFluidSpecWeight();
wLev = 10.; % level of the water table
domain.state.data.pressure = gamma_w*(z-wLev);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme.
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(simParam,domain);
% Solver = FCSolver(domain,'SaveRelError',true,'SaveBStepInf',true);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()