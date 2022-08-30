close all;
clear;
clc;

% --------------------------- READING MESH ---------------------------
% Getting parameters from the object mesh created with "Mesh.m"
mesh = Mesh();

% Starting a timer
tic;

% Setting the input file name
fileName = 'mesh.msh';

% Calling a function from the class "Mesh" that import mesh data
mesh.importGMSHmesh(fileName);

% Calling a function from the class "Mesh" for centroid calculation
mesh.finalize();

% Reading time
Reading_mesh = toc;
% ------------------------- END READING MESH -----------------------------

% ---------------------------- MESH DATA  --------------------------------
% Getting data from the object created with Mesh.m

% Element VTK type:
elemType = mesh.cellVTKType;

% Element material:
elemMAT = mesh.cellTag;

% Total number of element of the mesh
nElem = mesh.nCells;

% Total number of surfaces of the mesh
nFace = mesh.nSurfaces;

% Total number of nodes of the mesh
nTotNode = mesh.nNodes;

% Number of node for each element
nNode = mesh.cellNumVerts;

% Coordinates of each node
nodeCoords = mesh.coordinates;

% Topology of each element
elemTopol = mesh.cells;

% Topology of each surface
surfTopol = mesh.surfaces;

% Cells' centroids
cellCentroid = mesh.cellCentroid;
% -------------------------- END MESH DATA -------------------------------

%------------------------------- MATERIALS -------------------------------
% Starting a timer
tic;

% Setting the input file name
fileName = 'materials.dat';

% Creation of an object of "Materials"
mat = Materials(fileName);

% Calling the "getMaterial" function from the class "Materials" 
elas = mat.getMaterial('elas');

% Reading time
Generate_materials = toc;
%----------------------------- END MATERIALS ----------------------------

%-------------------------------- ELEMENTS ------------------------------
% Starting a timer
tic;

% Creation of an object of "Elements"
elems = Elements(nNode,nodeCoords,nElem,elemTopol,surfTopol,cellCentroid);
hexa = elems.getElement(elemType);

% Calling the function of the class that calculates the lenghts of sides of
% the element
hexa.getSide();

% Calling the function of the class that reorders the elements' topology
%hexa.reOrderTopol();

% Calling the function of the class that calculates the area of the faces 
% of the elements
%hexa.getArea();

% Calling the function of the class that calculates the volume of the
% element
%hexa.getVolume();

% Reading time
Generate_elements = toc;
%---------------------------- END ELEMENTS -------------------------------

%--------------------------- BOUNDARY CONDITIONS -------------------------
% Starting a timer
tic;

% DIRICHLET
% Setting input file
fileName = 'DirichletNode.dat';

%Creation of an object of "Boundaries" (general boundary condition)
bound = Boundaries(fileName);

% Defining boundary conditions in the input file as Dirichlet's boundary
% conditions for the nodes of the mesh (BCIdentifier = nodeDir)
nodedisp = bound.getBC('nodeDir');
% Creation of the object "nodedisp" (boundary conditions for node displacements)
nodedisp.NodeBoundary();

% NEUMANN
% Setting input file
fileName = 'NeumannNode.dat';

%Creation of an object of "Boundaries" (general boundary condition)
bound = Boundaries(fileName);

% Defining boundary conditions in the input file as Neumann's boundary
% conditions for the nodes of the mesh (BCIdentifier = nodeNeu)
nodeforce = bound.getBC('nodeNeu');
% Creation of the object "nodeforce" (boundary conditions for node forces)
nodeforce.NodeBoundary()

% Reading time
Generate_BC = toc;
%------------------------- END BOUNDARY CONDITIONS -----------------------

%------------------------------- ASSEMBLY --------------------------------
% Starting a timer
tic;

% Function to assembly global stiffness matrix
% If needed for constitutive law, specify stress value after "hexa"
K = assemblyKGauss(nTotNode, elemMAT, mat, hexa);
autoval = eigs(K);
%spy(K)
% Reading time
AssemblyK_time = toc;

% Function to assembly boundary conditions
[K,f] = assemblyBC(nTotNode, K, nodedisp, nodeforce);
% Reading time
AssemblyBC_time = toc;
%----------------------------- END ASSEMBLY ------------------------------

% Printing time 
fprintf('Time to read data from the mesh input file %.3f [s]\n', Reading_mesh);
fprintf('Time to generate object "Materials" %.3f [s]\n', Generate_materials);
fprintf('Time to generate object "Elements" %.3f [s]\n', Generate_elements);
fprintf('Time to generate object "Boundaries" %.3f [s]\n', Generate_BC);
fprintf('Time to assembly global stiffness matrix %.3f [s]\n', AssemblyK_time);
fprintf('Time to assembly boundary conditions %.3f [s]\n', AssemblyBC_time);

%--------------------------- SYSTEM SOLVING ------------------------------
% Starting a timer
tic;

u = K\f;
% Reading time
Solving_time = toc;
%------------------------- END SYSTEM SOLVING ----------------------------

% Printing time 
fprintf('Time to solve %.3f [s]\n', Solving_time);

%---------------------- DISPLACEMENT EVALUATION --------------------------
u_x = zeros(nTotNode,1);
u_y = zeros(nTotNode,1);
u_z = zeros(nTotNode,1);
for i = 1:nTotNode
  u_x(i) = u(i*3-2);
  u_y(i) = u(i*3-1);
  u_z(i) = u(i*3);
end

