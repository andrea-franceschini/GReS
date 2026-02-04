function mesh = structuredMesh(NX,NY,NZ,Lx,Ly,Lz)
% Return cartesian structured hexahedral grid

if nargin == 6
  % uniformly spaced grid
  x = linspace(Lx(1),Lx(2),NX+1);
  y = linspace(Ly(1),Ly(2),NY+1);
  z = linspace(Lz(1),Lz(2),NZ+1);
  nnx = NX + 1;
  nny = NY + 1;
  nnz = NZ + 1;
elseif nargin == 3
  % node spacing given as input input
  x = NX;
  y = NY;
  z = NZ;
  nnx = length(x);
  nny = length(y);
  nnz = length(z);
end



[X,Y,Z] = ndgrid(x,y,z);
coord = [X(:), Y(:), Z(:)];

% linear node numbering
nid = reshape(1:(nnx)*(nny)*(nnz), nnx, nny, nnz);

n1 = nid(1:end-1,1:end-1,1:end-1);
n2 = nid(2:end,1:end-1,1:end-1);
n3 = nid(1:end-1,2:end,1:end-1);
n4 = nid(2:end,2:end,1:end-1);
n5 = nid(1:end-1,1:end-1,2:end);
n6 = nid(2:end,1:end-1,2:end);
n7 = nid(1:end-1,2:end,2:end);
n8 = nid(2:end,2:end,2:end);

topol = [n1(:),n2(:),n4(:),n3(:),n5(:),n6(:),n8(:),n7(:)];

% boundary surfaces

% bottom: tag 1
n = nid(:,:,1);
fb = FaceGrid(n,nnx);

% top: tag 2
n = nid(:,:,end);
ft = FaceGrid(n,nnx);

% south: tag 3
n = nid(:,1,:);
fs = FaceGrid(n,nnx);

% north: tag 4
n = nid(:,end,:);
fn = FaceGrid(n,nnx);

% west
n = nid(1,:,:);
fw = FaceGrid(n,nny);

% east
n = nid(end,:,:);
fe = FaceGrid(n,nny);


surf = [fb; ft; fs; fn; fw; fe];

surfTag = [
    ones(size(fb,1),1) * 1
    ones(size(ft,1),1) * 2
    ones(size(fs,1),1) * 3
    ones(size(fn,1),1) * 4
    ones(size(fw,1),1) * 5
    ones(size(fe,1),1) * 6
];

mesh = Mesh();
mesh.coordinates = coord;
mesh.cells = topol;
mesh.surfaces = surf;
mesh.nCells = size(topol,1);
mesh.nSurfaces = size(surf,1);
mesh.surfaces = surf;
mesh.surfaceTag = surfTag;
mesh.cellTag = ones(mesh.nCells,1);
mesh.cellNumVerts = 8*ones(mesh.nCells,1);
mesh.surfaceNumVerts = 4*ones(mesh.nSurfaces,1);
mesh.cellVTKType = 12*ones(mesh.nCells,1);
mesh.surfaceVTKType = 9*ones(mesh.nSurfaces,1);
mesh.nSurfaceTag = 6;
mesh.nCellTag = 1;
mesh.nNodes = size(mesh.coordinates,1);
mesh.nDim = 3;

end


function f = FaceGrid(nodes,n)

f = [nodes(1:end-n-1)',...
     nodes(2:end-n)',...
     nodes(n+1:end-1)',...
     nodes(n+2:end)'];

f = f(:,[2 1 3 4]);
end

