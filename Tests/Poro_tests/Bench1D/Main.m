close all;
clear;
%
% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("Poromechanics_FEM");
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "simParam.dat";
simParam = SimulationParameters(fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'Bench1D_hexa2.msh';
%
% Import mesh data into the Mesh object
topology.importGMSHmesh(fileName);
topology.coordinates(:,3) = topology.coordinates(:,3)/5;
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'materialsList.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(fileName);
%
%------------------------------ ELEMENTS -----------------------------
%
% Define Gauss points
GaussPts = Gauss(12,2,3);
%
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);
%
% Create an object of the "Faces" class and process the face properties
faces = Faces(topology,model);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
%------------------------ BOUNDARY CONDITIONS ------------------------
% Setting input file
fileName = ["dirNodBotFace_hexa2.dat", "neuSurfTopFace_hexa2.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Some preprocessing stuff
pre = PreProc(grid,mat,GaussPts);
%
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat,pre,GaussPts);
%
% Create and set the print utility
printUtils = OutState(model,grid,'outTime.dat');
%
% ---------------------------- SOLUTION -------------------------------
%
% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam,grid,mat,pre,bound, ...
        printUtils,resState,GaussPts);
%
% Solve the problem
[simState] = NonLinearLoop(NSolv);
%
% Finalize the print utility
printUtils.finalize()
%
% Clean-up
delete(bound);