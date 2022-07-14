close all;
clear;
clc;

% Neumann and Dirichlet boundary conditions are defined in separeted input 
% files,however it is possibile to put together these two input files and generate
% different input files based on the type of BC (node, surface..) instead

%------------------------------  DIRICHLET  ------------------------------
% Setting input file
fileName = 'dirNode.dat';

% Creation of the object "Boundaries" (general boundary condition) 
bound = Boundaries(fileName);
% Starting a timer
tic;

% Defining boundary conditions in the input file as Dirichlet's boundary
% conditions for the nodes of the mesh (BCIdentifier = nodeDir)
nodedisp = bound.getBC('nodeDir');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
nodedisp.NodeBoundary()

% Reading the stopwatch timer
tD = toc;
%------------------------------  NEUMANN  -------------------------------
% Setting input file
fileName = 'neuNode.dat';

%Creation of an object of "Boundaries" (general boundary condition)
bound = Boundaries(fileName);
% Starting a timer
tic;

% Defining boundary conditions in the input file as Neumann's boundary
% conditions for the nodes of the mesh (BCIdentifier = nodeNeu)
nodeforce = bound.getBC('nodeNeu');
%  Creation of the object "nodeforce" (boundary conditions for node forces)
nodeforce.NodeBoundary()

% Reading the stopwatch timer
tN = toc;
% Printing t1 = time to read 
fprintf('Time to read and set Dirichlet BC %.3f [s]\n', tD);
fprintf('Time to read and set Neumann BC %.3f [s]\n', tN);

