close all;
clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

%% -------------------------- SET THE PHYSICS -------------------------
% model = ModelType("SinglePhaseFlow_FVTPFA");
model = ModelType("SinglePhaseFlow_FEM");

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = strcat(input_dir,'simParam.xml');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
fileName = strcat(input_dir,'Mesh/Fault.msh');

% Import mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%% ----------------------------- MATERIALS -----------------------------
% Set the input file name
fileName = strcat(input_dir,'materialsList.dat');
% fileName = strcat(input_dir,'materials.xml');

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%% ------------------------------ ELEMENTS -----------------------------
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model,topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

%% ----------------------- DOF Manager -----------------------------
% Degree of freedom manager
dofmanager = DoFManager(topology,model);

% Create and set the print utility
printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
    'folderName','Outputs');

%% ----------------------- Boundary Condition -----------------------------
% Creating and Appling boundaries conditions.
physics = "SinglePhaseFlow";
% physics = "VariablySaturatedFlow";
cond = struct('name',[],'type',[],'field',[],'values',[],'times',[]);
cond(1).name = 'BoundA';
cond(1).type = 'Dir';
cond(1).field = "latY0";
cond(1).times = 0.;
cond(1).values = 1e6;
cond(2).name = 'BoundB';
cond(2).type = 'Dir';
cond(2).field = "latYM";
cond(2).times = 0.;
cond(2).values = 1e5;

fileName = setBoundaryC('Inputs',grid,cond,physics);
% fileName = ["Inputs/BC_BoundA.dat","Inputs/BC_BoundB.dat"];
bound = Boundaries(fileName,model,grid);

%% ----------------------- Discretizer -----------------------------
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

%% ----------------------- Initial Condition -----------------------------
% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = 0;%1.e5;
% domain.state.data.potential(:) = domain.state.data.pressure+ mat.getFluid().getFluidSpecWeight()*topology.cellCentroid(:,3);

%% ----------------------- Solver -----------------------------
% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(domain,'SaveRelError',true,'SaveBStepInf',true);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()