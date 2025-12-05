close all;
clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

createBC = true;
% createBC = false;

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("SinglePhaseFlow_FVTPFA");

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
if createBC
  physics = "SinglePhaseFlow";
  cond = struct('name',[],'type',[],'field',[],'values',[],'times',[]);
  cond(1).name = 'BoundA';
  cond(1).type = 'Spg';
  cond(1).field = "latY0";
  cond(1).times = 0.;
  % cond(1).values = 1e6;
  cond(1).values = 3.5;
  % cond(1).times = [0. 1. 2. 3.];
  % cond(1).values = [1e6 1e6 1e6 1e6];

  % cond(2).name = 'Top';
  cond(2).name = 'BoundB';
  cond(2).type = 'Spg';
  cond(2).field = "latYM";
  cond(2).times = 0.;
  % cond(2).values = 1e5;
  cond(2).values = 3.;
  % cond(2).times = [0. 10. 20. 30.] ;
  % cond(2).values = [1e5 2e5 1e5 2e5];

  % cond(2).name = 'BoundB';
  % cond(2).type = 'Spg';
  % cond(2).field = "latYM";
  % cond(2).times = 0.;
  % cond(2).values = 1;

  cond(3).name = 'Top';
  cond(3).type = 'Neu';
  cond(3).field = "latZM";
  cond(3).times = 0.;
  cond(3).values = 0.;

  cond(4).name = 'Bot';
  cond(4).type = 'Neu';
  cond(4).field = "latZ0";
  cond(4).times = 0.;
  cond(4).values = 0.;

  % cond(5).name = 'LatA';
  % cond(5).type = 'Neu';
  % cond(5).field = "latXM";
  % cond(5).times = 0.;
  % cond(5).values = 0.;
  % 
  % cond(6).name = 'LatB';
  % cond(6).type = 'Neu';
  % cond(6).field = "latX0";
  % cond(6).times = 0.;
  % cond(6).values = 0.;

  fileName = setBoundaryC('Inputs',grid,cond,physics);
else
  fileName = ["Inputs/BC_BoundA.dat","Inputs/BC_BoundB.dat","Inputs/BC_Top.dat"];
end
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
z = elems.mesh.cellCentroid(:,3);
gamma_w = getFluid(mat).getFluidSpecWeight();
wLev = 1; % level of the water table
domain.state.data.pressure = gamma_w*(wLev-z);
% domain.state.data.pressure(:) = 1.e5;

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