clear
close all
% testing the local stiffness matrix assembly of an hexa and its partition
% in 6 tetra

% Using Gauss Points for exact integration, we investigate the influence of
% grid geometry on the stiffness matrix arising
% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("Poromechanics_FEM");
%
%
% ------------------------------  MESH -------------------------------
%
% Set the input file name
% BASE HEXA and TETRA meshes
hexamsh = Mesh();
hexamsh.importGMSHmesh('Hexa_1.msh');

tetramsh = Mesh();
tetramsh.importGMSHmesh('Tetra_1.msh');

% get fixed list of bounded nodes and loaded nodes


% list of nodes of top and bottom surfaces
bottom = find(hexamsh.coordinates(:,3) == 0);
top = find(hexamsh.coordinates(:,3) ~= 0);
dof = repelem(3*bottom',3);
dofBC = dof + repmat(-2:0,1,length(bottom));
% Uncomment the following two lines to obtain vertical linear benchmark
% dofNew = [top*3-1 top*3-2];
% dofBC = [dofBC (dofNew(:))'];

dofLoad = 3*top;

%----------------------------- MATERIALS -----------------------------
% Material properties
E = 1000000;
nu = 0.3;
D = zeros(6);
D([1 8 15]) = 1-nu;
D([2 3 7 9 13 14]) = nu;
D([22 29 36]) = (1-2*nu)/2;
D = (E/((1+nu)*(1-2*nu)))*D;

% Define Gauss points
GaussPts = Gauss(12,3,3);

dx_vals = [1 2 5 10 20 50 1000];
errNorm = zeros(length(dx_vals),1);
i = 0;
for dx = dx_vals
    %
    i = i+1;
    % Modify coordinates of mesh objects to deform the base geometry
    hexamsh.coordinates(hexamsh.coordinates(:,1) ~= 0 ,1) = dx;
    tetramsh.coordinates(tetramsh.coordinates(:,1) ~= 0 ,1) = dx;
    hexamsh.coordinates(hexamsh.coordinates(:,2) ~= 0 ,2) = dx;
    tetramsh.coordinates(tetramsh.coordinates(:,2) ~= 0 ,2) = dx;
    % Create an object of the "Elements" class and process the element properties
    elemHexa = Elements(hexamsh,GaussPts);
    elemTetra = Elements(tetramsh);
    facesHexa = Faces(model,hexamsh);
    facesTetra = Faces(model,tetramsh);
    % Create an object of the "Faces" class and process the face properties
    %
    % compute stiffness matrix
    K_hexa = computeStiffMat(hexamsh,elemHexa,D,GaussPts);
    K_tetra = computeStiffMat(tetramsh,elemTetra,D,GaussPts);
    F_hexa = zeros(24,1);
    F_tetra = zeros(24,1);
    aInfHexa = computeAreaNod(hexamsh,facesHexa);
    aInfTetra = computeAreaNod(tetramsh,facesTetra);
    F_hexa(dofLoad) = -ones(length(dofLoad),1).*aInfHexa(top);
    F_tetra(dofLoad) = -ones(length(dofLoad),1).*aInfTetra(top);
    fDir = zeros(length(dofBC),1);
    [K_hexa,F_hexa] = applyDir(dofBC, fDir, K_hexa, F_hexa);
    [K_tetra,F_tetra] = applyDir(dofBC, fDir, K_tetra, F_tetra);
    u_hexa = K_hexa\F_hexa;
    u_tetra = K_tetra\F_tetra;
    rel = (u_hexa-u_tetra)./u_hexa;
    rel(isnan(rel)) = 0; rel(isinf(rel))= 0;
    errNorm(i) = norm(rel,2);
end
% % %%
% diff = cell(4,1);
% relNorm = zeros(4,1);
% for i = 1:4
%     diff{i} = u{i}(u{i}~=0)-u{i+3}(u{i+3}~=0);
%     relNorm(i) = norm(diff{i}./u{i}(u{i}~=0));
% end



% %%
% K_hexa = {K{1}; K{2}; K{3}};
% K_tetra = {K{4}; K{5}; K{6}};
% relDiff1 = abs(K_hexa{1} - K_tetra{1})./K_hexa{1};
% relDiff2 = abs(K_hexa{2} - K_tetra{2})./K_hexa{2};
% relDiff3 = abs(K_hexa{3} - K_tetra{3})./K_hexa{3};
% n1 = norm(relDiff1,"fro");
% n2 = norm(relDiff2,"fro");
% n3 = norm(relDiff3,"fro");


% Loading the cube with a vertical load and checking the relative difference
% of the displacement solution


%fac = (E*(1-nu))/((1+nu)*(1-2*nu));
fac = D(end);
anal = 1/fac;
