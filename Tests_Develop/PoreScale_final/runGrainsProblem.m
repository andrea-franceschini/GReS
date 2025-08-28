isModelReady = true; % flag to avoid reprocessing domains and interfaces

if ~isModelReady

  fprintf('\n Processing grain grid \n')

  %read the grain mesh ( only volumes)
  grainMesh = Mesh();
  grainMesh.importMesh('Mesh/many_grains.vtk');
  grainMesh.cellTag = grainMesh.cellTag + 1;
  grainMesh.nCellTag = grainMesh.nCellTag + 1;

  mult_type = 'P0';

  % generate surface mesh for each grain
  grainMesh = setSurfaces(grainMesh,mult_type);
  elems = Elements(grainMesh,4);

  % get grain connectivity
  % produce a pair of connectivity, each represents a master/slave pair to
  % introduce mortar interfaces
  conn = getGrainConnectivity(grainMesh);

  % ensure that if a grain is not connected to any other grain, then at least
  % lies on the dirichlet boundary (unit cube bounding box)
  % just need to check that at least one node has one coordinate equal to max
  % or min coordinates of the domain
  isFree = find(~ismember(1:grainMesh.nSurfaceTag,unique(conn)));

  for i = isFree
    c = grainMesh.coordinates(grainMesh.surfaces(grainMesh.surfaceTag==i,:),:);
    if isempty(c)
      continue
    end
    tol = 1e-3;

    isDir =  any(abs(c-min(c))<tol,"all") & any(abs(c-max(c))<tol,"all");
    assert(isDir,'Grain %i is floating \n',i)
  end




  % write interface file
  interfFile = 'Domains/interfaces.xml';
  writeInterfaceFile(interfFile,conn,mult_type);


  simParam = SimulationParameters('simParam.dat');

  % setup model
  model = ModelType("Poromechanics_FEM");

  faces = Faces(model,grainMesh);

  grid = struct('topology',grainMesh,'cells',elems,'faces',faces);

  material = Materials(model,'Materials/MaterialsList.dat');

  dof = DoFManager(grainMesh,model);


  printUtils = OutState(model,grainMesh,'outTime.dat','writeVtk',true,'folderName','OUT_p0');

  linSyst = Discretizer(model,simParam,dof,grid,material);



  % creare domains struct (BCs initially empty)
  domains = struct( ...
    'DomainName',         'Grains', ...
    'ModelType',          model, ...
    'Grid',               grid, ...
    'Material',           material, ...
    'DoFManager',         dof, ...
    'BoundaryConditions', [], ...
    'OutState',           printUtils, ...
    'Discretizer',        linSyst ...
    );

  fprintf('\n Processing mortar interfaces \n')

  % setup mortar interfaces
  [interfaces,domains] = Mortar.buildInterfaceStruct(interfFile,domains);

  %% Process boundary conditions
  processBC_grains(grainMesh,interfaces,press_grain);
  domains.BoundaryConditions = Boundaries(["BCs/pressureX.dat","BCs/pressureY.dat","BCs/pressureZ.dat",...
    "BCs/fixZ.dat","BCs/fixY.dat","BCs/fixX.dat","BCs/fixAll.dat"],domains.ModelType,domains.Grid);


  % save mat files
  save('Interface_grain2grain.mat',"interfaces");
  save('domains_grain2grain.mat',"domains");

else

  % interfaces and domains already processed

  interf = load('Interface_grain2grain.mat');
  interfaces = interf.interfaces;

  dom = load('domains_grain2grain.mat');
  domains = dom.domains;

  for i = 1:numel(interfaces)
    interfaces{i}.solvers(1) = domains.Discretizer;
    interfaces{i}.solvers(2) = domains.Discretizer;
  end


end
%%
fprintf('\n Running simulation \n')

simParam = SimulationParameters('simParam.dat');

solver = MultidomainFCSolver(simParam,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();


































































% load grain interfaces from mat files


% 
% for i = 1:numel(interfaces_g2g)
%   interfaces_g2g{i}.solvers(1) = domains.Discretizer;
%   interfaces_g2g{i}.solvers(2) = domains.Discretizer;
% end
% 
% grainMesh = domains.Grid.topology;
% 
% nodesInterf = [];
% % get list of surfaces of grain2grain interface
% for i = 1:numel(interfaces_g2g)
%   nodesInterf = [nodesInterf;(interfaces_g2g{i}.mesh.local2glob{1})';(interfaces_g2g{i}.mesh.local2glob{2})'];
% end
% 
% nodesInterf = unique(nodesInterf);
% 
% isContactSurf = ~any(~ismember(grainMesh.surfaces,nodesInterf),2);
% 
% % get list of surfaces to apply fluid pressure:
% % exclude boundary surfaces and grain2grain surfaces
% 
% tol = 1e-3;
% % get surfaces having only boundary nodes
% isBoundSurf = any([abs(grainMesh.surfaceCentroid-0.15)<tol,abs(grainMesh.surfaceCentroid-0.85)<tol],2);
% 
% nodesTop = find(abs(grainMesh.coordinates(:,3)-0.85)<tol);
% nodesBot = find(abs(grainMesh.coordinates(:,3)-0.15)<tol);
% nodesNorth = find(abs(grainMesh.coordinates(:,2)-0.85)<tol);
% nodesSouth = find(abs(grainMesh.coordinates(:,2)-0.15)<tol);
% nodesEast = find(abs(grainMesh.coordinates(:,1)-0.15)<tol);
% nodesWest = find(abs(grainMesh.coordinates(:,1)-0.85)<tol);
% nodesBound = unique([nodesWest;nodesEast;nodesSouth;nodesNorth;...
%   nodesBot;nodesTop]);
% 
% 
% 
% isFluidSurf = ~any([isBoundSurf isContactSurf],2); 
% % get id of nodes on the fluid surfaces
% nodeFluid = unique(grainMesh.surfaces(isFluidSurf,:));
% nodeFluid = nodeFluid(~ismember(nodeFluid,nodesBound));
% 
% % get nodal forces (pressure*area*averagenodenormal)
% nodeNormal = zeros(grainMesh.nNodes,3);
% nodeArea = zeros(grainMesh.nNodes,1);
% for i = reshape(find(isFluidSurf),1,[])
%   nid = grainMesh.surfaces(i,:);
%   nodeArea(nid) = nodeArea(nid)+grainMesh.surfaceArea(i)/3;
%   coords = grainMesh.coordinates(nid,:);
%   v1 = coords(1,:)-coords(2,:);
%   v2 = coords(1,:)-coords(3,:);
%   n = cross(v1,v2)/norm(cross(v1,v2));
%   nodeNormal(nid,:) = nodeArea(nid).*n; 
% end
% 
% nodeNormal = nodeNormal(nodeFluid,:);   % already considering the node area!
% nodeArea = nodeArea(nodeFluid);
% 
% press_nodes = press_grain(nodeFluid).*nodeNormal;
% pX = press_nodes(:,1);
% pY = press_nodes(:,2);
% pZ = press_nodes(:,3);
% 
% 
% 
% surfTags = [];
% % spot boundary grains lying on the boundary (here we must fix all nodes)
% % fix all associated dofs 
% % get all effective surfaceTag (included in the input of the interfaces)
% for i = 1:numel(interfaces_g2g)
%   surftaglist = [interfaces_g2g{i}.mesh.msh(1).surfaceTag;...
%     interfaces_g2g{i}.mesh.msh(2).surfaceTag];
%   surfTags = [surfTags;unique(surftaglist)];
% end
% 
% surfTags = unique(surfTags);
% floatTags = ~ismember(1:grainMesh.nSurfaceTag,surfTags);
% 
% nodeFloating = unique(grainMesh.surfaces(ismember(grainMesh.surfaceTag,[find(floatTags) 138]),:));
% 
% % modify bcs of domain structure
% 
% 
% writeBCfiles('BCs/fixZ','NodeBC','Dir',{'Poromechanics','z'},'FixedZ',0,0,[nodesTop;nodesBot]);
% writeBCfiles('BCs/fixX','NodeBC','Dir',{'Poromechanics','x'},'FixedX',0,0,[nodesWest;nodesEast]);
% writeBCfiles('BCs/fixY','NodeBC','Dir',{'Poromechanics','y'},'FixedY',0,0,[nodesSouth;nodesNorth]);
% writeBCfiles('BCs/fixAll','NodeBC','Dir',{'Poromechanics','x','y','z'},'FixedAll',0,0,nodeFloating);
% writeBCfiles('BCs/pressureX','NodeBC','Neu',{'Poromechanics','x'},'LoadX',0,0,nodeFluid,-pX);
% writeBCfiles('BCs/pressureY','NodeBC','Neu',{'Poromechanics','y'},'LoadY',0,0,nodeFluid,-pY);
% writeBCfiles('BCs/pressureZ','NodeBC','Neu',{'Poromechanics','z'},'LoadZ',0,0,nodeFluid,-pZ);
% 
% domains.BoundaryConditions = Boundaries(["BCs/pressureX.dat","BCs/pressureY.dat","BCs/pressureZ.dat",...
%   "BCs/fixZ.dat","BCs/fixY.dat","BCs/fixX.dat","BCs/fixAll.dat"],domains.ModelType,domains.Grid);
% 
% end
% 
% 
% 
% simParam = SimulationParameters('simParam.dat');
% 
% 
% solver = MultidomainFCSolver(simParam,domains,interfaces);
% solver.NonLinearLoop();
% solver.finalizeOutput();


%% AUXILIARY ROUTINES

function msh = setSurfaces(msh,mult_type)
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
  % get unique faces (lie on the boundary of the grain)
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

  msh.surfaces = [msh.surfaces;...
    surfs];
  msh.surfaceTag = [msh.surfaceTag;...
    i*ones(N,1)];
end


M = max(msh.coordinates,[],'all');
m = min(msh.coordinates,[],'all');

% centr = zeros(size(msh.surfaces,1),3);
% 
% for jj = 1:size(centr,1)
%   xc = mean(msh.coordinates(msh.surfaces(jj,:),1));
%   yc = mean(msh.coordinates(msh.surfaces(jj,:),2));
%   zc = mean(msh.coordinates(msh.surfaces(jj,:),3));
%   centr(jj,:) = [xc yc zc];
% end
% tol = 1e-3;
% id = any([any(abs(centr-M)<tol,2), any(abs(centr-m)<tol,2)],2);
% % finally remove boundary surfaces
% msh.surfaces(id,:) = [];
% msh.surfaceTag(id) = [];


% better safe then sorry: check that no surfaces left have all three nodes in
% the boundary
tol = 1e-3;
list = false(size(msh.surfaces,1),1);
for i=1:size(msh.surfaces,1)
  n = msh.surfaces(i,:);
  c = msh.coordinates(n,:);
  isSurfBound = all(any([any(abs(c-M)<tol,2), any(abs(c-m)<tol,2)],2));
  if isSurfBound
    list(i) = true;
  end
end
msh.surfaces(list,:) = [];
msh.surfaceTag(list) = [];


if ~strcmp(mult_type,'P0')

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


function conn = getGrainConnectivity(msh,varargin)
conn = [];
ng = msh.nSurfaceTag;
if isempty(varargin)
    list = 1:ng; 
  else
    list = varargin{1};
  end
  for i = list
    fprintf('Processing connectivity of grain %i \n',i)
    for j = i+1:ng
      % rough screen using bounding box expaned
      if ~doBBIntersect(msh,i,j)
        continue
      elseif checkIntersect(msh,i,j)
        % proper contact detection
        conn = [conn;...
                i,j];
      end
    end
  end
end

function out = doBBIntersect(msh,i1,i2)
  % check intersection of two grains by comparing the bounding boxes
  % (expanded by 5% of their size for proximity contact)

  c = msh.coordinates(msh.surfaces(msh.surfaceTag == i1,:),:);
  min1 = [min(c(:,1)) min(c(:,2)) min(c(:,3))];
  max1 = [max(c(:,1)) max(c(:,2)) max(c(:,3))];
  d = abs(max1-min1);
  min1 = min1 - 0.025*d;
  max1 = max1 + 0.025*d;
  c = msh.coordinates(msh.surfaces(msh.surfaceTag == i2,:),:);
  min2 = [min(c(:,1)) min(c(:,2)) min(c(:,3))];
  max2 = [max(c(:,1)) max(c(:,2)) max(c(:,3))];
  d = abs(max2-min2);
  min2 = min2 - 0.025*d;
  max2 = max2 + 0.025*d;

  out = all(max1 >= min2) & all(max2 >= min1);

end

function out = checkIntersect(msh,i1,i2)
 out = false;
 s1 = getSurfaceMesh(msh,i1);
 s2 = getSurfaceMesh(msh,i2);
 if any([isempty(s1) isempty(s2)])
   return
 end
 cs = ContactSearching(s1,s2);
 out = sum(cs.elemConnectivity,'all')>0;
end

function writeInterfaceFile(fname,conn,mult_type)

slaves = unique(conn(:,1));
N = numel(slaves);
% setup base structure
interfBase.Type = 'MeshTying';
interfBase.Quadrature.typeAttribute = 'SegmentBased';
interfBase.Quadrature.nGPAttribute = 7;
interfBase.Quadrature.nIntAttribute = 5;
if strcmp(mult_type,'P0')
  interfBase.Stabilization.typeAttribute = "Jump";
end
interfBase.Multiplier.typeAttribute = mult_type;
interfBase.Physics = 'Poromechanics';
interfBase.Master.idAttribute = num2str(1);
interfBase.Slave.idAttribute = num2str(1);

Interface = repmat(interfBase,N,1);

for i = 1:N
  % get list of masters for each slave
  id = conn(:,1)==slaves(i);
  m = num2str(reshape(conn(id,2),1,[]));
  Interface(i).Master.surfaceTagAttribute = m;
  Interface(i).Slave.surfaceTagAttribute = num2str(slaves(i));
end

str.Interface = Interface;

writestruct(str,fname);
end

