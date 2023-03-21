close all;
clear;
addpath('../../Code/read');
addpath('../../Code/boundaries');
addpath('../../Code/discretizer');
addpath('../../Code/elements');
addpath('../../Code/materials');
addpath('../../Code/gauss');
addpath('../../Code/preproc');
addpath('../../Code/state');
addpath('../../Code/write');
addpath('../../Code/nonlinearsolver');
addpath('../../Code/modeltype');
addpath('../../Code/outstate');
addpath('../../Code/simparam');

% -------------------------- SET THE PHYSICS -------------------------
%
model = ModelType('Poromechanics');
%
% ------------------- READING SIMULATION PARAMETERS ------------------
%
fileName = "SimParam.dat";
simParam = SimulationParameters(fileName);
%
% --------------------------- READING MESH ---------------------------
% Getting parameters from the object mesh created with "Mesh.m"
mesh = Mesh();
%
% Setting the input file name
fileName = 'mesh.msh';
%
% Calling a function from the class "Mesh" that import mesh data
mesh.importGMSHmesh(fileName);
%
[mesh.cells] = reOrderTopol(mesh.cells,mesh.coordinates,mesh.nCells);
%
%------------------------------- MATERIALS -------------------------------
%
% Setting the input file name
fileName = 'materials.dat';

% Creation of an object of "Materials"
mat = Materials(fileName);
%
%-------------------------------- ELEMENTS ------------------------------
% Define Gauss points
GaussPts = Gauss(12,2,3);
%
% Creation of an object of "Elements"
elems = Elements(mesh,GaussPts);
hexa = elems.getElement(mesh.cellVTKType);
%
%--------------------------- BOUNDARY CONDITIONS -------------------------
%
% Setting input file
fileName = ["dirNode.dat", "neuNode.dat"];

%Creation of an object of "Boundaries" (general boundary condition)
bound = Boundaries(fileName);
%
%------------------------------- ASSEMBLY --------------------------------
%
% Computing some pre-processing stuff
pre = PreProc(mesh,mat,GaussPts);
% Define the initial reservoir state
resState = State(model,mesh,hexa,mat,pre,GaussPts);
%
% Setting up the print utility 
printUtils = OutState(model,mesh,'outTime.dat');
%
% Define the name of the BC that are to be imposed
BCName = ["nodeDir", "nodeNeu"];
%
% Nonlinear solver in time
% Setup
NSolv = NonLinearSolver(model,simParam,mesh,hexa,mat,pre,bound,BCName, ...
        printUtils,resState,GaussPts);
% Performing the loop
[simState] = NonLinearLoop(NSolv);
%
% Clean-up
delete(bound);