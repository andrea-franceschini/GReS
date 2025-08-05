%wrap fluid model 
% model = ModelType("SinglePhaseFlow_FEM");
% dofFluid = DoFManager(gridFluid.topology,model);
% material = Materials(model,'Materials/materialsFlow.dat');
% linSystFluid = Discretizer(model,simParam,dofFluid,gridFluid,material);
% 
% domains(1) = struct( ...
%   'DomainName',         'Fluid', ...
%   'ModelType',          model, ...
%   'Grid',               gridFluid, ...
%   'DoFManager',         dofFluid, ...
%   'Discretizer',        linSystFluid ...
%   );
% 
% 
% % wrap grains model
% grainMesh = Mesh();
% grainMesh.importMesh('Mesh/many_grains.vtk');
% 
% % merge all cell tags into unique one
% grainMesh.cellTag = ones(grainMesh.nCells,1);
% grainMesh = setSurfacesGrain(grainMesh,'P0');
% model = ModelType("SinglePhaseFlow_FEM");
% 
% elems = Elements(grainMesh,4);
% faces = Faces(model,grainMesh);
% grid = struct('topology',grainMesh,'cells',elems,'faces',faces);
% dof = DoFManager(grainMesh,model);
% material = Materials(model,'Materials/materialsFlow.dat');
% linSyst = Discretizer(model,simParam,dof,grid,material);
% domains(2) = struct( ...
%   'DomainName',         'Grains', ...
%   'ModelType',          ModelType("SinglePhaseFlow_FEM"), ...
%   'Grid',               grid, ...
%   'DoFManager',         dof, ...
%   'Discretizer',        linSyst ...
%   );
% 
% % setup mortar interfaces
% [interfaces,domains] = Mortar.buildInterfaceStruct('Domains/fluid2grain.xml',domains);

%%

interf = load("Interface.mat","interfaces");
interfaces{1} = interf.interfaces{1};

%interfaces{1}.computeMortarMatrices();

interfaces{1}.computeMat();

E = interfaces{1}.E;

% get active nodes in the flow domain (remove nodes belonging to floating
% cells)
actEnts = interfaces{1}.solvers(1).dofm.fields.isEntActive;


% interpolate cell pressure into node using node area
mshFlow = interfaces{1}.solvers(1).grid.topology;
mshGrain = interfaces{1}.solvers(2).grid.topology;
cellFlow =interfaces{1}.solvers(1).grid.cells;
cellGrain = interfaces{1}.solvers(2).grid.cells;
press_node = zeros(mshFlow.nNodes,1);
vol_node = zeros(mshFlow.nNodes,1);
for i = 1:mshFlow.nCells
  nId = mshFlow.cells(i,:);
  volNod = findNodeVolume(cellFlow,i);
  press_node(nId) = press_node(nId) + press(i)*volNod; 
  vol_node(nId) = vol_node(nId) + volNod; 
end


area_node = zeros(mshGrain.nNodes);
for i = 1:mshGrain.nSurfaces
  nId = mshGrain.surfaces(i,:);
  areaNod = findNodeArea(cellGrain,i);
  area_node(nId) = area_node(nId) + areaNod'; 
end



press_node = press_node./vol_node;

% interpolate pressure from flow mesh to grain surfaces
press_grain = E*press_node(actEnts);

% get id of elements where D matrix has crazy value (due to diagonal
% approximation for partially integrated slave)
id = any([press_grain > 11 , press_grain < -1],2); 

idD = find(ismember(interfaces{1}.mesh.local2glob{2},find(id)));

a = area_node(id);
E(id,:) = interfaces{1}.Jmaster{1}(idD,:).*a;
press_grain = E*press_node(actEnts);


% modify the grain mesh removing surfaces that lies on the external
% boundary % print and processing purposes


tol = 1e-3;
% get surfaces having only boundary nodes
boundSurf = any([abs(mshGrain.surfaceCentroid-0.15)<tol,abs(mshGrain.surfaceCentroid-0.85)<tol],2);
mshGrain.surfaceTag(boundSurf) = 2;
mshGrain.nSurfaceTag = mshGrain.nSurfaceTag + 1;
[mshGrain,gnodes] = mshGrain.getSurfaceMesh(1);

press_grain = press_grain(unique(gnodes));
  















function msh = setSurfacesGrain(msh,mult_type)
% set surfaces with proper tag for each grain
for i = 1:msh.nCellTag
  % extract grain i
  c = msh.cells(msh.cellTag==i,:);
  % define faces for each grain
  s = [c(:,[1 2 3]);...
    c(:,[1 3 4]);...
    c(:,[1 2 4]);...
    c(:,[2 3 4])];
  [s,i1] = sort(s,2);
  % get unique faces (lie on the boundary)
  [s_u,i2,i3] = unique(s,'rows');
  isFaceBound = accumarray(i3,1) == 1;
  id = i1(i2,:);
  id = id(isFaceBound,:);
  N = sum(isFaceBound);
  [k1,~] = find(id' == 1);
  [k2,~] = find(id' == 2);
  [k3,~] = find(id' == 3);
  s1 = sub2ind([N 3],(1:N)',k1);
  s2 = sub2ind([N 3],(1:N)',k2);
  s3 = sub2ind([N 3],(1:N)',k3);
  s_u = s_u(isFaceBound,:);
  surfs = [s_u(s1) s_u(s2) s_u(s3)];

    % remove surfaces lying on the bounding box (avoid issues with bcs)
    % get centroid of each cell
    msh.surfaces = [msh.surfaces;...
      surfs];
    msh.surfaceTag = [msh.surfaceTag;...
      i*ones(N,1)];
end

if ~strcmp(mult_type,'P0')

  centr = zeros(size(msh.surfaces,1),3);

  for jj = 1:size(centr,1)
    xc = mean(msh.coordinates(msh.surfaces(jj,:),1));
    yc = mean(msh.coordinates(msh.surfaces(jj,:),2));
    zc = mean(msh.coordinates(msh.surfaces(jj,:),3));
    centr(jj,:) = [xc yc zc];
  end
  tol = 1e-3;
  M = max(msh.coordinates,[],'all');
  m = min(msh.coordinates,[],'all');
  id = any([any(abs(centr-M)<tol,2), any(abs(centr-m)<tol,2)],2);
  msh.surfaces(id,:) = [];
  msh.surfaceTag(id) = [];

  % further remove any surface containing one node on the top boundary
  z = zeros(size(msh.surfaces,1),3);
  for kk = 1:size(msh.surfaces,1)
    z(kk,:) = msh.coordinates(msh.surfaces(kk,:),3);
  end

  id = any(abs(z-M)<tol,2);
  msh.surfaces(id,:) = [];
  msh.surfaceTag(id) = [];
end

msh.nSurfaces = numel(msh.surfaceTag);
msh.surfaceVTKType = 5*ones(msh.nSurfaces,1);
msh.nSurfaceTag = max(msh.surfaceTag);
msh.surfaceNumVerts = 3*ones(msh.nSurfaces,1);
end
