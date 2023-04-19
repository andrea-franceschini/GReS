close all;
clear;
%
% -------------------------- SET THE PHYSICS -------------------------
%
model = ModelType('Poromechanics');
%
% --------------------- READING SIMULATION SETTINGS ------------------
%
fileName = "simParam.dat";
simParam = SimulationParameters(fileName);
%
% --------------------------- READING MESH ---------------------------
%
% Creation of an object with class "Mesh"
mesh = Mesh();
%
% Setting the input file name
fileName = 'mesh.msh';
%
% Calling a function from the class "Mesh" that imports mesh data
mesh.importGMSHmesh(fileName);
%
%------------------------------- MATERIALS -------------------------------
%
% Setting the input file name
fileName = 'materials.dat';
%
% Creation of an object of "Materials"
mat = Materials(fileName);
%
%-------------------------------- ELEMENTS ------------------------------
%
% Creation of an object of "Elements"
elems = Elements(mesh);
tetra = elems.getElement(mesh.cellVTKType);
%
%--------------------------- BOUNDARY CONDITIONS -------------------------
%
% Setting input file
fileName = ["dirNode.dat","neuNode.dat"];

%Creation of an object of "Boundaries" (general boundary condition)
bound = Boundaries(fileName);
%
%------------------------------- ASSEMBLY --------------------------------
%
% Computing some pre-processing stuff
pre = PreProc(mesh,mat);
%
% Define the initial reservoir state
resState = State(model,mesh,tetra,mat,pre);
%
% Setting up the print utility 
printUtils = OutState(model,mesh,'outTime.dat');
%
% Define the name of the BC that are to be imposed
BCName = ["nodeDir","nodeNeu"];
%
% Nonlinear solver in time
% Setup
NSolv = NonLinearSolver(model,simParam,mesh,tetra,mat,pre,bound,BCName, ...
        printUtils,resState);
% Performing the loop
[simState] = NonLinearLoop(NSolv);
%
% Finalize the print utility
printUtils.finalize()
%
% Clean-up
delete(bound);
