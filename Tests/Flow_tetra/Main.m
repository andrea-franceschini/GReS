close all;
clear;

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("SinglePhaseFlow_FEM");
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "SimParam.dat";
simParam = SimulationParameters(fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'mesh_25x100x5.msh';
%
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'materials.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(fileName);
%
%------------------------------ ELEMENTS -----------------------------
%
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',[]);
%
%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = "dirNode_2wells.dat";
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Some preprocessing stuff
pre = PreProc(grid,mat);
%
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat,pre);
%
% Create and set the print utility
printUtils = OutState(model,grid,'outTime.dat');
%
% Print the reservoir initial state
printUtils.printState(resState);
%
% ---------------------------- SOLUTION -------------------------------
%
NSolv = NonLinearSolver(model,simParam,grid,mat,pre,bound,printUtils,resState);
%
% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%
delete(bound);