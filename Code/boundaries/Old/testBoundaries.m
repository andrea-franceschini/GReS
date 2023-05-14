close all;
clear;

% Neumann and Dirichlet boundary conditions are defined in separeted input 
% files, however it is possibile to put together these two input files and generate
% different input files based on the type of BC (node, surface..) instead

%------------------------------  DIRICHLET  ------------------------------
% Setting input file
fileName = ["dirNodePoro.dat","neuNodePoro.dat","dirNodeFlow.dat", ...
            "neuNodeFlow.dat","volForceFlow.dat","volForcePoro.dat"];
% fileName = ["dirNodePoro.dat","neuNodePoro.dat","dirNodeFlow.dat"];

% Creation of the object "Boundaries" (general boundary condition) 
bound = Boundaries(fileName);
% Starting a timer
tic;

% Defining boundary conditions in the input file as Dirichlet's boundary
% conditions for the nodes of the mesh (BCIdentifier = nodeDir)
nodeDirP = bound.getBC('nodeDirP');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
nodeDirP.boundPhysics
nodeDirP.boundType
nodeDirP.boundDof
nodeDirP.boundVal
nodeDirP.timeInt
%
nodeDirF = bound.getBC('nodeDirF');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
nodeDirF.boundPhysics
nodeDirF.boundType
nodeDirF.boundDof
nodeDirF.boundVal
nodeDirF.timeInt
%
nodeNeuP = bound.getBC('nodeNeuP');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
nodeNeuP.boundPhysics
nodeNeuP.boundType
nodeNeuP.boundDof
nodeNeuP.boundVal
nodeNeuP.timeInt
%
nodeNeuF = bound.getBC('nodeNeuF');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
nodeNeuF.boundPhysics
nodeNeuF.boundType
nodeNeuF.boundDof
nodeNeuF.boundVal
nodeNeuF.timeInt
%
volForceF = bound.getBC('volForceF');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
volForceF.boundPhysics
volForceF.boundDof
volForceF.boundVal
volForceF.timeInt
%
volForceP = bound.getBC('volForceP');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
volForceP.boundPhysics
volForceP.boundDof
volForceP.boundVal
volForceP.timeInt

% Reading the stopwatch timer
tD = toc;
BCName = ["nodeDirP", "nodeDirF", "nodeNeuP", "nodeNeuF", "volForceF", "volForceP"];
for i=1:length(BCName)
cond = getBC(bound,BCName(i));
t = 3.5;
Boundaries.checkAndUpdateCond(cond,t);
cond.boundPhysics
% cond.boundType
cond.boundDof
cond.boundVal
cond.timeInt
end

%------------------------------  NEUMANN  -------------------------------
% Setting input file
% fileName = 'neuNode.dat';

%Creation of an object of "Boundaries" (general boundary condition)
% bound = Boundaries(fileName,'Neu');
% Starting a timer
tic;

% Defining boundary conditions in the input file as Neumann's boundary
% conditions for the nodes of the mesh (BCIdentifier = nodeNeu)
% nodeNeu = bound.getBC('nodeNeu');
% %  Creation of the object "nodeforce" (boundary conditions for node forces)
% nodeNeu.boundPhysics
% nodeNeu.boundType
% nodeNeu.boundDof
% nodeNeu.boundVal

% Reading the stopwatch timer
tN = toc;
% Printing t1 = time to read 
fprintf('Time to read and set Dirichlet BC %.3f [s]\n', tD);
fprintf('Time to read and set Neumann BC %.3f [s]\n', tN);
% delete(bound);
% clear bound

