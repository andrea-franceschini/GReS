clear
close all
clc


scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);
%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('Faulted aquifer model \n')
fprintf('___________________\n\n')

% get mesh file
mshStr = load("cornerPoint.mat","mesh");
mesh = mshStr.mesh;
mesh.nCellTag = 1;

% set surfaces and return location of wells
wellsId = setAquiferSurfaces(mesh);



% Set physical models 
model = ModelType("Poromechanics_FEM");

% Set parameters of the simulation
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);


% Create an object of the Materials class and read the materials file
fileName = 'materialsList.dat';
mat = Materials(model,fileName);


% Create an object of the "Elements" class and process the element properties
gaussOrder = 2;
elems = Elements(mesh,gaussOrder);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, mesh);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',mesh,'cells',elems,'faces',faces);
%
dofmanager = DoFManager(mesh, model);

% Create and set the print utility
printUtils = OutState(model,mesh,'outTime.dat','folderName','Faulted_aquifer','flagMatFile',false);


% testing surfaces
% surf3 = getSurfaceMesh(mesh,3);
% surf4 = getSurfaceMesh(mesh,4);
% surf5 = getSurfaceMesh(mesh,5);
% surf6 = getSurfaceMesh(mesh,6);
% surf7 = getSurfaceMesh(mesh,7);
% surf8 = getSurfaceMesh(mesh,8);
% 
% plotFunction(surf3,"surf3Test",ones(surf3.nNodes,1));
% plotFunction(surf4,"surf4Test",ones(surf4.nNodes,1));
% plotFunction(surf5,"surf5Test",ones(surf5.nNodes,1));
% plotFunction(surf6,"surf6Test",ones(surf6.nNodes,1));
% plotFunction(surf7,"surf7Test",ones(surf7.nNodes,1));
% plotFunction(surf8,"surf8Test",ones(surf8.nNodes,1));

% Create an object of the "Boundaries" class 
% write BC files
setBCfiles(mesh,wellsId);
fileName = ["BCs/fix_X.dat", "BCs/fix_Y.dat", "BCs/fix_Z.dat", "BCs/source.dat"];
bound = Boundaries(fileName,model,grid);


% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

interfFile = 'Domains/interface.xml';

% set verbosity 
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domain);


interfaces{1}.solvers(2).simparams.setVerbosity(2);

solver = MultidomainFCSolver(domain,interfaces);
%solver = ActiveSetContactSolver(domain,interfaces,5);

%solver.simParameters.setBackstepSkipFlag(1);

solver.NonLinearLoop();
solver.finalizeOutput();

%% plot profiles of multipliers along vertical axis (avoid opening paraview)
%plotStep(solver.results,2);


function wellsId = setAquiferSurfaces(mesh)

% find external surfaces on the grid and fix them

% find all faces of the grid
% find external faces
% set surfaces with proper tag for each grain

c = mesh.cells;


% full face list (unsorted, keeps identity)
s_full = [c(:,[1 2 3 4]); ...
          c(:,[1 2 6 5]); ...
          c(:,[2 3 7 6]); ...
          c(:,[3 4 8 7]); ...
          c(:,[1 4 8 5]); ...
          c(:,[5 6 7 8])];

idFace = repelem((1:6)', mesh.nCells, 1);

% for uniqueness, sort each row
[s_sorted, ~] = sort(s_full, 2);
[~, i2, i3] = unique(s_sorted, 'rows');

% boundary mask
isFaceBound = accumarray(i3,1) == 1;

% recover boundary faces and their ids
s_bnd     = s_full(i2(isFaceBound), :);   % proper vertex ordering
idFace_bnd = idFace(i2(isFaceBound));     % which face (1..6)

% discard surfaces already present (like faults)
id = all(ismember(s_bnd,mesh.surfaces),2);

s_bnd = s_bnd(~id,:);
idFace_bnd = idFace_bnd(~id);



% add surface tags
mesh.surfaces = [mesh.surfaces;...
  s_bnd];
mesh.surfaceTag = [mesh.surfaceTag;...
  idFace_bnd+2];



% finalize the mesh object
mesh.nSurfaces = numel(mesh.surfaceTag);
mesh.surfaceNumVerts = 4*ones(mesh.nSurfaces,1);
mesh.surfaceVTKType = 9*ones(mesh.nSurfaces,1);

% locate wells id 
% map location of cells in i-j-k index to cell id
i1 = [10,10,7];
i2 = [10,30,7];

wellsId = zeros(2,1);
wellsId(1) = sub2ind([34,40,14], i1(1), i1(2), i1(3));
wellsId(2) = sub2ind([34,40,14], i2(1), i2(2), i2(3));






end





% 

