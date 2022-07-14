close all;
clear;
clc;

% --------------------------- READING MESH ---------------------------
% % Getting parameters from the object mesh created with "Mesh.m"
% mesh = Mesh();
% 
% % Starting a timer
% tic;
% 
% % Setting the input file name
% fileName = 'mesh.msh';
% 
% % Calling a function from the class "Mesh" that import mesh data
% mesh.importGMSHmesh(fileName);
% 
% % Reading time
% Reading_mesh = toc;
% ------------------------- END READING MESH -----------------------------
%
% ---------------------------- MESH DATA  --------------------------------
% Getting data from the object created with Mesh.m
%
% Element type:
% elemType = mesh.cellVTKType;
%
% Total number of element of the mesh
% nElem = mesh.nCells;
% 
% Total number of surfaces of the mesh
% nFace = mesh.nSurfaces;
% 
% Total number of nodes of the mesh
% nNode = mesh.nNodes
% 
% Coordinates of each node
% nodeCoords = mesh.coordinates;
%
% Topology of each element
% elemTopol = mesh.cells;
%
% Topology of each surface
% surfTopol = mesh.surfaces;
% -------------------------- END MESH DATA -------------------------------

% DATI PER IL SOLO TEST DI FUNZIONAMENTO DELLA CLASSE:
elemType = 10;
nNode = 6;
nElem = 2;
elemTopol = [1 3 2 5;3 4 2 6];
surfTopol = [3 5 1;
             3 2 1;
             3 2 5;
             1 2 5;
             4 6 3;
             3 6 2;
             4 2 6;
             3 4 2];
nodeCoords = [1 1 0;
              2 3 0;
              3 1 0;
              4 3 0;
              2 2 2;
              3 2.5 1];

% ----------------------------- ELEMENTS ---------------------------------
% Starting a timer
tic;

% Creation of an object of "Elements"
elem = Elements(nNode,nodeCoords,nElem,elemTopol,surfTopol);          
tetra = elem.getElement(elemType);

% Calling the function of the class that calculates the coefficient for the
% stiffness matrix
tetra.getCoefficient();

% Calling the function of the class that calculates the volume of the
% element
tetra.getVolume()

% Calling the function of the class that calculates the matrix of the 
% derivates of the elements
tetra.getDerivates();

% Calling the function of the class that calculates the area of the faces 
% of the elements
tetra.getArea()
% ---------------------------- END ELEMENTS ------------------------------

% Reading the stopwatch timer
tE = toc;
% Printing time
%fprintf('Time to read data from input file %.3f [s]\n', Reading_mesh);
fprintf('Time to generate object "element" %.3f [s]\n', tE);
