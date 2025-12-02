function [mesh, pressures] = runFlowSimulation(mesh,gridDims,nRock)

fprintf('Faulted aquifer model - flow simulation \n')
fprintf('___________________\n\n')

wellsId = setAquiferMesh(mesh,gridDims,nRock);

fprintf('ID of cells with wells: %i    %i \n',wellsId(1),wellsId(2))

% Set physical models 
model = ModelType("SinglePhaseFlow_FVTPFA");

% Set parameters of the simulation
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);

% Create an object of the Materials class and read the materials file
fileName = 'Materials/materialsListFlow.dat';
mat = Materials(model,fileName);

% Create an object of the "Elements" class and process the element properties
gaussOrder = 1;
elems = Elements(mesh,gaussOrder);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, mesh);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',mesh,'cells',elems,'faces',faces);
%
dofmanager = DoFManager(mesh, model);

% Create and set the print utility
printUtils = OutState(model,mesh,'outTime.dat','folderName','Faulted_aquifer_flow','flagMatFile',true);


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
setBCflow(mesh,wellsId);
fileName = ["BCs/fix_press.dat", "BCs/wells_press.dat"];
bound = Boundaries(fileName,model,grid);


% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

solver = FCSolver(domain);
solver.NonLinearLoop();
domain.outstate.finalize()

pressures = printUtils.results.expPress;

end




function setBCflow(mesh,wellsId)

writeBCfiles('BCs/fix_press','SurfBC','Dir','SinglePhaseFlow','fix_press',0,0,mesh,[5,7,8]);
writeBCfiles('BCs/wells_press','ElementBC','Dir','SinglePhaseFlow','wells',[0,3,7,10],[0,120,440,750], wellsId);

end


function wellsId = setAquiferMesh(mesh,gridDims,nRock)

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

i1 = [nRock+round(0.6*(gridDims(1)-nRock)),round(0.45*gridDims(2)),round(0.5*gridDims(3))];
i2 = [nRock+round(0.6*(gridDims(1)-nRock)),round(0.55*gridDims(2)),round(0.5*gridDims(3))];

wellsId = zeros(2,1);
wellsId(1) = sub2ind(gridDims, i1(1), i1(2), i1(3));
wellsId(2) = sub2ind(gridDims, i2(1), i2(2), i2(3));

% assign cell tags for rock, clay sand
rock = [1,nRock; 1,gridDims(2); 1,gridDims(3)];
clay = [nRock+1,gridDims(1); 1,gridDims(2); round(0.6*gridDims(3))+1,gridDims(3)];
sand = [nRock+1,gridDims(1); 1,gridDims(2); 1,round(0.6*gridDims(3))];

cells = cell(3,1);

ic = 0;

for mat = {rock,clay,sand}
  m = mat{1};
  cellList = zeros(prod(m(:,2)-m(:,1)),1);
  kk = 0;
  for i=m(3,1):m(3,2)
    for j=m(2,1):m(2,2)
      for k=m(1,1):m(1,2)
        kk = kk+1;
        cellList(kk) = sub2ind(gridDims, k,j,i); 
      end
    end
  end
  cells{ic+1} = cellList;
  ic = ic+1;
end

for i = 1:3
  mesh.cellTag(cells{i}) = i;
end

mesh.nCellTag = max(mesh.cellTag);

end
