close all;
input_dir = 'Inputs/';

typeFlow = "Saturated";
%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create the Mesh object
topology = Mesh();

% Import mesh data into the Mesh object
topology.importMesh(fullfile(input_dir,'Mesh','dam.msh'));

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Create an object of the "Faces" class and process the face properties
faces = Faces(topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Creating boundaries conditions.
bcs = ["boundariesDir.xml","boundariesSPG.xml"];
bound = Boundaries(fullfile(input_dir,bcs(2)),grid);

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility for the solution
printUtils = OutState(fullfile(input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound);

switch typeFlow
  case "Saturated"
    domain.addPhysicsSolver(fullfile(input_dir,'solverSaturated2.xml'));
  case "Unsaturated"
    domain.addPhysicsSolver(fullfile(input_dir,'solverUnsaturated.xml'));
end

% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = 1.e5;
% % % % domain.state.data.pressure(:) = -10*9.8066e3;
% % % % domain.state.data.pressure(:) = 1e5;

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme.
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function

% Solve the problem
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();









% % % 
% % % % cond(1).name = 'wLev';
% % % % cond(1).type = 'Dir';
% % % % cond(1).field = "topA";
% % % % cond(1).times = 0.;
% % % % % cond(1).values = -0.75*9.8066e3;
% % % % cond(1).values = 1e3;
% % % 
% % % % cond(1).name = 'wLev';
% % % % cond(1).type = 'Spg';
% % % % cond(1).field = "topA";
% % % % cond(1).times = 0.;
% % % % % cond(1).values = -0.75*9.8066e3;
% % % % cond(1).values = 23.5;
% % % 
% % % % cond(2).name = 'Ux0';
% % % % cond(2).type = 'Dir';
% % % % cond(2).field = "latX0";
% % % % cond(2).times = 0.;
% % % % cond(2).values = -10*9.8066e3;
% % % % cond(2).name = 'Ux0';
% % % % cond(2).type = 'Spg';
% % % % cond(2).field = "latX0";
% % % % cond(2).times = 0.;
% % % % cond(2).values = 20.;
% % % 
% % % % cond(3).name = 'Uz0';
% % % % cond(3).type = 'Neu';
% % % % cond(3).field = "bot";
% % % % cond(3).times = 0.;
% % % % cond(3).values = 0;
% % % % 
% % % % cond(3).name = 'NeuB';
% % % % cond(3).type = 'Neu';
% % % % cond(3).field = "topB";
% % % % cond(3).times = 0.;
% % % % cond(3).values = 0.;
% % % % 
% % % % cond(4).name = 'NeuC';
% % % % cond(4).type = 'Neu';
% % % % cond(4).field = "topC";
% % % % cond(4).times = 0.;
% % % % cond(4).values = 0.;
% % % 
% % % 
% % % % cond(1).name = 'Ux0';
% % % % cond(1).type = 'Dir';
% % % % cond(1).field = "latX0";
% % % % cond(1).times = 0.;
% % % % cond(1).values = 1e6;
% % % 
% % % % cond(2).name = 'UxM';
% % % % cond(2).type = 'Dir';
% % % % cond(2).field = "latXM";
% % % % cond(2).times = 0.;
% % % % cond(2).values = 1e5;
% % % 
% % % 
% % % % cond(1).name = 'Ux0';
% % % % cond(1).type = 'Spg';
% % % % cond(1).field = "latX0";
% % % % cond(1).times = 0.;
% % % % cond(1).values = 20.;
% % % 
% % % cond(1).name = 'wLev';
% % % cond(1).type = 'Spg';
% % % cond(1).field = "topA";
% % % cond(1).times = 0.;
% % % cond(1).values = 4.5;
% % % 
% % % 