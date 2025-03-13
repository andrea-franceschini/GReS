close all;
clear;

scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);

% simple mechanical model 
% A cube of size 1x1x1m is fixed in the bottom face and a load is applied
% at the top face. Both a load in vertical and horizontal direction is
% considered. Note that Drucker-Prager plasticity will not activate if only
% a vertical load is considered.
% Load magnitude is just indicative and must be tuned to correctly activate
% the non linear model.

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("Poromechanics_FEM");
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
%fileName = 'Mesh/cube.msh';
fileName = 'Mesh/column.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
%
fileName = 'Materials/MaterialsList.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);
%
%------------------------------ ELEMENTS -----------------------------
%
GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% create a DoFManager object
dofmanager = DoFManager(topology, model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

% Build a structure storing variable fields at each time step
state = linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','vtkOutput');

%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["BCs/bottom_fix.dat", "BCs/load_z.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);

% Set initial stress state
nu = mat.db(1).ConstLaw.nu;
z0 = -1000; % average depth
gamma = 12000; % N/m^3
M1 = 0.50; %nu / (1-nu);
M2 = 0.50; %nu / (1-nu);
cells = dofmanager.getFieldCells("Poromechanics");
for el = cells'
  switch topology.cellVTKType(el)
    case 10
      z = z0 + elems.cellCentroid(el,2);
      sigmaz = gamma * z;
      sigmax = M1 * sigmaz;
      sigmay = M2 * sigmaz;
      state.iniStress(el,1) = sigmax;
      state.iniStress(el,2) = sigmay;
      state.iniStress(el,3) = sigmaz;
  end
end
state.conv.stress = state.iniStress;

%
% ---------------------------- SOLUTION -------------------------------
%
% Create the object handling the (nonlinear) solution of the problem
solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);
[simState, endState] = solver.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%
