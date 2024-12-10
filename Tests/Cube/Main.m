close all;
clear;

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
simParam = SimulationParameters(model,fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'Mesh/Cube.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);


% UTILITY TO WRITE BC FILES (check Code/boundaries folder in github feature/moretto)
% writeBCfiles('bottom_fix','NodeBC','Dir','Poro',["x","y","z"],'bottom_fix',0,0,topology,1);
% writeBCfiles('load_z','SurfBC','Neu','Poro',"z",'top_vertLoad',[0 1],[0 -100],topology,2);
%writeBCfiles('load_x','SurfBC','Neu','Poro',"x",'top_shearLoad',[0 1],[0 10],topology,2);
%%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
%
%fileName = 'Materials/MaterialsList.dat';
fileName = 'Materials/materialsList2.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);
%
%Terzaghi_analytical(topology, mat, 10)
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
%
% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["BCs/bottom_fix.dat","BCs/load_x.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid,dofmanager);

% Create a State object that store the solution vectors during solution
% loop
resState = State(model,grid,mat,GaussPts);

% Create and set the print utility
printUtils = OutState(model,mat,grid,'outTime.dat');
%
% Print model initial state
printUtils.printState(resState);
%
% ---------------------------- SOLUTION -------------------------------
%
% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,resState,GaussPts);
%
% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%