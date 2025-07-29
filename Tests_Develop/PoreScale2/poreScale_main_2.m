clear
close all
clc

%% Grain mesh pre-processing
% only 1 domain grid is define, with different physical volumes (cellTags)
% interacting

% a preprocessing routine identifies the surface of each grain, and assign
% a tag to each one. 

% read the grain mesh ( only volumes)
grainMesh = Mesh();
grainMesh.importMesh('Mesh/many_grains.vtk');
grainMesh.cellTag = grainMesh.cellTag + 1;
grainMesh.nCellTag = grainMesh.nCellTag + 1;


% generate surface mesh for each grain
grainMesh = setSurfaces(grainMesh);
elems = Elements(grainMesh,4);

% get grain connectivity
% produce a pair of connectivity, each represents a master/slave pair to
% introduce mortar interfaces



interfFile = 'Domains/interfaces.xml';

plotFunction(grainMesh,'domain',ones(grainMesh.nNodes,1))

% for i = 1:grainMesh.nSurfaceTag
%   mshTmp = getSurfaceMesh(grainMesh,i);
%   plotFunction(mshTmp,strcat('out_',num2str(i)),ones(mshTmp.nNodes,1))
% end

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

% create xml file for connecting domains
writeInterfaceFile(interfFile,conn);


%%

nodesTop = find(abs(grainMesh.coordinates(:,3)-0.85)<tol);
nodesBot = find(abs(grainMesh.coordinates(:,3)-0.15)<tol);
nodesNorth = find(abs(grainMesh.coordinates(:,2)-0.85)<tol);
nodesSouth = find(abs(grainMesh.coordinates(:,2)-0.15)<tol);
nodesEast = find(abs(grainMesh.coordinates(:,1)-0.15)<tol);
nodesWest = find(abs(grainMesh.coordinates(:,1)-0.85)<tol);


writeBCfiles('BCs/moveTop','NodeBC','Dir',{'Poromechanics','z'},'MovingTop',0,-0.01,nodesTop);
writeBCfiles('BCs/fixBot','NodeBC','Dir',{'Poromechanics','z'},'FixedBottom',0,0,nodesBot);
writeBCfiles('BCs/fixX','NodeBC','Dir',{'Poromechanics','x'},'FixedX',0,0,[nodesWest;nodesEast]);
writeBCfiles('BCs/fixY','NodeBC','Dir',{'Poromechanics','y'},'FixedY',0,0,[nodesSouth;nodesNorth]);

%% MODEL
simParam = SimulationParameters('simParam.dat');

% setup model
model = ModelType("Poromechanics_FEM");

faces = Faces(model,grainMesh);

grid = struct('topology',grainMesh,'cells',elems,'faces',faces);

material = Materials(model,'Materials/MaterialsList.dat');

dof = DoFManager(grainMesh,model);

bc = Boundaries(["BCs/moveTop.dat","BCs/fixBot.dat", "BCs/fixX.dat", "BCs/fixY.dat"],model,grid);

printUtils = OutState(model,grainMesh,'outTime.dat','writeVtk',true,'folderName','OUT_p0');

linSyst = Discretizer(model,simParam,dof,grid,material);


domains = struct( ...
  'DomainName',         'Grains', ...
  'ModelType',          model, ...
  'Grid',               grid, ...
  'Material',           material, ...
  'DoFManager',         dof, ...
  'BoundaryConditions', bc, ...
  'OutState',           printUtils, ...
  'Discretizer',        linSyst ...
  );

% setup mortar interfaces
[interfaces,domains] = Mortar.buildInterfaceStruct(interfFile,domains);

solver = MultidomainFCSolver(simParam,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();

















%% AUXILIARY ROUTINES

function msh = setSurfaces(msh)
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

function writeInterfaceFile(fname,conn)

slaves = unique(conn(:,1));
N = numel(slaves);
% setup base structure
interfBase.Type = 'MeshTying';
interfBase.Quadrature.typeAttribute = 'SegmentBased';
interfBase.Quadrature.nGPAttribute = 7;
interfBase.Quadrature.nIntAttribute = 5;
%interfBase.Stabilization.typeAttribute = "Jump";
interfBase.Multiplier.typeAttribute = 'dual';
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


