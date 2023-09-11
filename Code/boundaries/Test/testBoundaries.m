close all;
clear;

addpath('../');

% Neumann and Dirichlet boundary conditions are defined in separeted input
% files, however it is possibile to put together these two input files and generate
% different input files based on the type of BC (node, surface..) instead

%------------------------------  DIRICHLET  ------------------------------
% Setting input file
fileName = ["dirNodePoro.dat","neuNodePoro.dat","dirNodeFlow.dat", ...
  "neuNodeFlow.dat","volForceFlow.dat","volForcePoro.dat"];

% Creation of the object "Boundaries" (general boundary condition)
bound = Boundaries(fileName);
% Starting a timer
tic;

% % Defining boundary conditions in the input file as Dirichlet's boundary
% % conditions for the nodes of the mesh (BCIdentifier = nodeDir)
%nodeDirP = bound.getBCData('nodeDirP');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
%nodeDirP

%
%nodeDirF = bound.getBCData('nodeDirF');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
%nodeDirF

%
%nodeNeuP = bound.getBCData('nodeNeuP');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
%nodeNeuP
%
%nodeNeuF = bound.getBCData('nodeNeuF');
%  Creation of the object "nodedisp" (boundary conditions for node displacements)
%nodeNeuF

t = 15;
% Reading the stopwatch timer
tD = toc;
BCName = ["nodeDirP", "nodeDirF", "nodeNeuP", "nodeNeuF", "volForceF", ...
  "volForceP"];
for i=1:length(BCName)
  vals = getVals(bound,BCName(i), t);
  dofs = getDofs(bound,BCName(i));
  [dofs,vals]
end

% Printing t1 = time to read
fprintf('Time to read and set BC %.3f [s]\n', tD);

bound.delete();
