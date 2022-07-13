close all;
clear;
clc;

% --------------------------- READING MESH ---------------------------
% % Getting parameters from the object mesh created with "Mesh.m"
% mesh = Mesh();
% 
% Starting a timer
  tic;
% 
% % Setting the input file name
% fileName = 'mesh.msh';
% 
% % Calling a function from the class "Mesh" that import mesh data
% mesh.importGMSHmesh(fileName);

%DATI DEGLI ELEMENTI DEL DOMINIO CREATO CON STRAUS, SE USASSI LA CLASSE
%MESH NON SAREBBE NECESSARIO ESPLICARE LE PROPERTIES (elemType ecc.)

%test assemblaggio matrice di rigidezza con più tetraedri
% Tipo di elemento VTK
elemType = 10;
% Numero di nodi di un elemento
nNode = 4;
% Numero totale di nodi nella mesh
nTotNode = 9;
% Numero totale di elementi nella mesh
nElem = 12;
% Matrice contenente la topologia degli elementi
elemTopol = [1	2	3	9
             3	4	1	9
             5	1	4	9
             4	8	5	9
             5	6	2	9
             2	1	5	9
             2	6	3	9
             7	3	6	9
             7	8	3	9
             4	3	8	9
             6	5	7	9
             8	7	5	9];
         
%surfTopol scorretta ma tanto per ora non serve
surfTopol = [3 1 4];

%Matrice contenente le coordinate dei nodi della mesh
nodeCoords = [4.00E+00	1.00E+01	0.00E+00
              4.00E+00	7.00E+00	0.00E+00
              7.00E+00	7.00E+00	0.00E+00
              7.00E+00	1.00E+01	0.00E+00
              4.00E+00	1.00E+01	5.00E+00
              4.00E+00	7.00E+00	5.00E+00
              7.00E+00	7.00E+00	5.00E+00
              7.00E+00	1.00E+01	5.00E+00
              5.50E+00	8.50E+00	2.50E+00];
%-------------------------- END READING MESH -----------------------------

%------------------------ READING MATERIAL -------------------------------
% Setting the input file name
fileName = 'materials.dat';

% Creation of an object of "Materials"
mat = Materials(fileName);

% Calling the "getMaterial" function from the class "Materials" 
elas = mat.getMaterial('elas');

% Calling the function of the class "Elastic" that calculates stiffness
% matrix for elastic materials
D = elas.getStiffnessMatrix();
%------------------------- END READING MATERIALS -------------------------

%------------------------ READING ELEMENTS -------------------------------
% Creation of an object of "Elements"
% elems = Elements(mesh);
elems = Elements(nNode,nodeCoords,nElem,elemTopol,surfTopol);
tetra = elems.getElement(elemType);

% Calling the function of the class that calculates the coefficient for the
% stiffness matrix
tetra.getCoefficient();

% Calling the function of the class that calculates the volume of the
% element
tetra.getVolume();

% Calling the function of the class that calculates the matrix of the 
% derivates for the single element
tetra.getDerivates();
%------------------------- END READING ELEMENTS --------------------------

%----------------------- READING BOUNDARY CONDITIONS ---------------------
% probabilmente si può fare un file di input unico
% DIRICHLET
% Setting input file
fileName = 'DirichletNode.dat';

%Creation of an object of "Boundaries"
bound = Boundaries(fileName);

nodedisp = bound.getBC('nodeDir');
nodedisp.NodeBoundary()

% NEUMANN
% Setting input file
fileName = 'NeumannNode.dat';

%Creation of an object of "Boundaries"
bound = Boundaries(fileName);

nodeforce = bound.getBC('nodeNeu');
nodeforce.NodeBoundary()
%------------------------- END BOUNDARY CONDITIONS -----------------------
Reading_time = toc;

%------------------------------- ASSEMBLY --------------------------------
% Starting a timer
tic;

% Function for assembly global stiffness matrix in three steps:
% 1) Assembly local stiffness matrices
% 2) Assembly boundary conditions
% 3) Assembly known term

[K,f] = assembly(nTotNode, D, tetra, nodedisp, nodeforce);
autoval = eigs(K);
spy(K)
%----------------------------- END ASSEMBLY ------------------------------

Assembly_time = toc;
% Printing time 
fprintf('Time to read %.3f [s]\n', Reading_time);
fprintf('Time to assembly %.3f [s]\n', Assembly_time);

%--------------------------- SYSTEM SOLVING ------------------------------
% Starting a timer
tic;

u = K\f;
f = K*u;
%------------------------- END SYSTEM SOLVING ----------------------------

Solving_time = toc;
% Printing time 
fprintf('Time to solve %.3f [s]\n', Solving_time);

