close all;
clear;
%
% -------------------------- SET THE PHYSICS -------------------------
%
model = ModelType('SinglePhaseFlow');
%
% ------------------- READING SIMULATION PARAMETERS ------------------
%
fileName = "simParam.dat";
simParam = SimulationParameters(fileName);
%
% --------------------------- READING MESH ---------------------------
%
% Getting parameters from the object mesh created with "Mesh.m"
mesh = Mesh();

% Setting the input file name
fileName = 'mesh_25x100x5.msh';

% Calling a function from the class "Mesh" that import mesh data
mesh.importGMSHmesh(fileName);
%
%------------------------------- MATERIALS -------------------------------
%
% Setting the input file name
fileName = 'materials.dat';

% Creation of an object of "Materials"
mat = Materials(fileName);
%
%-------------------------------- ELEMENTS ------------------------------
%
% Define Gauss points only if required, i.e., there is at least one
% hexahedral element
%
% Creation of an object of "Elements"
elems = Elements(mesh);
tetra = elems.getElement(mesh.cellVTKType);
%
%--------------------------- BOUNDARY CONDITIONS -------------------------
%
% Setting input file
% fileName = "dirNode_2wells.dat";
fileName = ["dirNode_1well.dat", "neuNode.dat"];

%Creation of an object of "Boundaries" (general boundary condition)
bound = Boundaries(fileName);
%
%------------------------------- ASSEMBLY --------------------------------
%
% Computing some pre-processing stuff
pre = PreProc(mesh,mat);
% Define the initial reservoir state
resState = State(model,mesh,tetra,mat,pre);
%
% Setting up the print utility 
printUtils = OutState(model,mesh,'outTime.dat');
printUtils.printState(resState);
%
% Define the name of the BC that are to be imposed
% BCName = "nodeDir";
BCName = ["nodeDir","nodeNeu"];
%
% Nonlinear solver in time
% Setup
NSolv = NonLinearSolver(model,simParam,mesh,tetra,mat,pre,bound,BCName, ...
        printUtils,resState);
% Performing the loop
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%
% Clean-up
delete(bound);
