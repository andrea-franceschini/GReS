close all;
clear;
clc;

% --------------------------- READING MESH ---------------------------
% Creation of an object with class "Mesh"
mesh = Mesh();

% Starting a timer
tic;

% Setting the input file name
fileName = 'mesh.msh';

% Calling a function from the class "Mesh" that imports mesh data
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


% Linear elastic isotropic material (matIdentifier = elas)
mat = Materials(fileName);
elas = mat.getMaterial('elas');

% Reading time
Generate_materials = toc;
%----------------------------- END MATERIALS ----------------------------

%-------------------------------- ELEMENTS ------------------------------
% Starting a timer
tic;

% Tetrahedrons 
elems = Elements(nNode,nodeCoords,nElem,elemTopol,surfTopol);
tetra = elems.getElement(elemType);

% Basis functions calculation
tetra.getBasisF();

% Calling the function of the class that calculates the volume of the
% element
tetra.getVolume();

% Basis function derivatives calculation 
tetra.getDerivatives();

% Reading time
Generate_elements = toc;
%---------------------------- END ELEMENTS -------------------------------

%--------------------------- BOUNDARY CONDITIONS -------------------------
% Starting a timer
tic;

% DIRICHLET
% Setting input file
fileName = 'dirNode.dat';

%Creation of an object of "Boundaries" (general boundary condition)
bound = Boundaries(fileName);

% Nodal Dirichlet boundary conditions (BCIdentifier = nodeDir)
nodedisp = bound.getBC('nodeDir');
% nodedisp = nodal displacements
nodedisp.NodeBoundary()

% NEUMANN
% Setting input file
fileName = 'neuNode.dat';

bound = Boundaries(fileName);

% Nodal Neumann boundary conditions (BCIdentifier = nodeNeu)
nodeforce = bound.getBC('nodeNeu');
% nodeforce = nodal forces
nodeforce.NodeBoundary()

% Reading time
Generate_BC = toc;
%------------------------- END BOUNDARY CONDITIONS -----------------------

%------------------------------- ASSEMBLY --------------------------------
% Starting a timer
tic;

% Function to assembly global stiffness matrix
K = assemblyK(nTotNode, elemMAT, mat, tetra);
autoval = eigs(K);

% Reading time
AssemblyK_time = toc;

% Function to assign boundary conditions
[K,f] = imposeBC(nTotNode, K, nodedisp, nodeforce);
% Reading time
ImposeBC_time = toc;
%----------------------------- END ASSEMBLY ------------------------------

% Printing time 
fprintf('Time to read data from the mesh input file %.3f [s]\n', Reading_mesh);
fprintf('Time to generate object "Materials" %.3f [s]\n', Generate_materials);
fprintf('Time to generate object "Elements" %.3f [s]\n', Generate_elements);
fprintf('Time to generate object "Boundaries" %.3f [s]\n', Generate_BC);
fprintf('Time to assembly global stiffness matrix %.3f [s]\n', AssemblyK_time);
fprintf('Time to impose boundary conditions %.3f [s]\n', ImposeBC_time);

%------------------------- LINEAR SYSTEM SOLVER --------------------------
% Starting a timer
tic;

u = K\f;
% Reading time
Solving_time = toc;
%------------------------- END SYSTEM SOLVER -----------------------------

%---------------------- DISPLACEMENT EVALUATION --------------------------
u_x = zeros(nTotNode,1);
u_y = zeros(nTotNode,1);
u_z = zeros(nTotNode,1);
for i = 1:nTotNode
  u_x(i) = u(i*3-2);
  u_y(i) = u(i*3-1);
  u_z(i) = u(i*3);
end

% ----------------------------- VTK RESULTS -----------------------------
addpath('../../write');

pointData3D = repmat(struct('name', 1, 'data', 1), 3, 1);
pointData3D(1).name = 'ux';
pointData3D(1).data = u_x;
pointData3D(2).name = 'uy';
pointData3D(2).data = u_y;
pointData3D(3).name = 'uz';
pointData3D(3).data = u_z;
cellData3D = repmat(struct('name', 1, 'data', 1), 0, 1);

V = VTKOutput(mesh);
V.writeVTKFile(0.0, pointData3D, [], [], []);

% ---------------------- ANALITYCAL SOLUTION ---------------------------
% Surface pressure
sigma_z = -0.2;  %MPa

% Vertical compressibility:
D = elas.getStiffnessMatrix();
K_z = D(3,3);

% Model geometry:
H = 1500;    %mm
z_coord = nodeCoords(:,3);

% Analitycal solution
u_real = zeros(nTotNode,1);
for i = 1:nTotNode
    u_real(i) = (sigma_z/K_z)*z_coord(i);
end

% ------------------------- 2-NORM ERROR ------------------------------
norm1 = norm(u_real - u_z);
norm2 = norm(u_real);
error = norm1/norm2;