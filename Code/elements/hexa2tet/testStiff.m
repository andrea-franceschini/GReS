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
fNames = ["Hexa_1.msh";
    "Hexa_2.msh";
    "Hexa_4.msh";
    "Hexa_5.msh";
    "Tetra_1.msh";
    "Tetra_2.msh";
    "Tetra_4.msh";
    "Tetra_5.msh"];

% get fixed list of bounded nodes and loaded nodes


F = zeros(24,1);
u = cell(length(fNames),1);

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

K = cell(length(fNames),1);

for i = 1:length(fNames)
    %
    % Import the mesh data into the Mesh object
    % Create the Mesh object
    topology = Mesh();
    f = convertStringsToChars(fNames(i));
    topology.importGMSHmesh(f);
    bottomBC = find(topology.coordinates(:,3) == 0);
    topLoad = find(topology.coordinates(:,3) ~= 0);
    % bottom surface fixed
    dof = repelem(3*bottomBC',3);
    dofBC = dof + repmat(-2:0,1,length(bottomBC));
    % lateral surface fixed in normal direction
    % dofNew = [topLoad*3-1 topLoad*3-2];
    % dofBC = [dofBC (dofNew(:))'];
    % free dofs (for error computation)
    freeDof = find(~ismember(1:24,dofBC));
    % loaded dofs
    dofLoad = 3*topLoad;
    %
    % Create an object of the "Elements" class and process the element properties
    elems = Elements(topology,GaussPts);
    faces = Faces(model,topology);
    % Create an object of the "Faces" class and process the face properties
    %
    % compute stiffness matrix
    K = computeStiffMat(topology,elems,D,GaussPts);
    F = zeros(24,1);
    aInf = computeAreaNod(topology,faces);
    F(dofLoad) = -ones(length(dofLoad),1).*aInf(topLoad);
    fDir = zeros(length(dofBC),1);
    [K,F] = applyDir(dofBC, fDir, K, F);
    u{i} = K\F;
end
% %%
diff = cell(4,1);
relNorm = zeros(4,1);
for i = 1:4
    diff{i} = u{i}(u{i}~=0)-u{i+3}(u{i+3}~=0);
    relNorm(i) = norm(diff{i}./u{i}(u{i}~=0));
end



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


fac = (E*(1-nu))/((1+nu)*(1-2*nu));
anal = 1/fac;
