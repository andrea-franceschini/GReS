function grid = structuredMesh(NX,NY,NZ,Lx,Ly,Lz)
% STRUCTUREDMESH  Build a Cartesian structured hexahedral mesh.
%
%   MESH = STRUCTUREDMESH(NX, NY, NZ, Lx, Ly, Lz) creates a uniform mesh
%   with NX x NY x NZ cells over the axis-aligned box defined by the
%   coordinate ranges Lx, Ly, Lz.
%
%   MESH = STRUCTUREDMESH(X, Y, Z) creates a mesh whose node positions are
%   given explicitly by the coordinate vectors X, Y, Z.
%
% -------------------------------------------------------------------------
% INPUTS (uniform spacing)
%   NX, NY, NZ  - number of cells along x, y, z  [positive integers]
%   Lx          - x extent: [x_min, x_max]        [1x2 double]
%   Ly          - y extent: [y_min, y_max]         [1x2 double]
%   Lz          - z extent: [z_min, z_max]         [1x2 double]
%
% INPUTS (prescribed spacing)
%   NX, NY, NZ  - node coordinate vectors along x, y, z  [double arrays]
%
% -------------------------------------------------------------------------
% OUTPUT
%   grid  - grid object with processed geometry
%
% -------------------------------------------------------------------------
% EXAMPLES
%   % Uniform 10x5x3 mesh on [0,1] x [0,0.5] x [0,0.3]
%   mesh = structuredMesh(10, 5, 3, [0 1], [0 0.5], [0 0.3]);

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
fb = FaceGrid(n);

% top: tag 2
n = nid(:,:,end);
ft = FaceGrid(n);


% south: tag 3
n = nid(:,1,:);
fs = FaceGrid(n);

% north: tag 4
n = nid(:,end,:);
fn = FaceGrid(n);

% west
n = nid(1,:,:);
fw = FaceGrid(n);

% east
n = nid(end,:,:);
fe = FaceGrid(n);


surf = [fb; ft; fs; fn; fw; fe];

surfTag = [
    ones(size(fb,1),1) * 1
    ones(size(ft,1),1) * 2
    ones(size(fs,1),1) * 3
    ones(size(fn,1),1) * 4
    ones(size(fw,1),1) * 5
    ones(size(fe,1),1) * 6
];

grid = Grid();
grid.nDim = 3;
grid.coordinates = coord;

nC = size(topol,1);
grid.cells.connectivity = topol;
grid.cells.VTKType = 12*ones(nC,1);
grid.cells.tag = ones(nC,1);
grid.cells.numVerts = 8*ones(nC,1);


nS = size(surf,1);
grid.surfaces.connectivity = surf;
grid.surfaces.VTKType = 9*ones(nS,1);
grid.surfaces.tag = surfTag;
grid.surfaces.numVerts = 4*ones(nS,1);

% finalize geometry
grid.setStructured;
processGeometry(grid);

end


function f = FaceGrid(nodes)

sz = size(nodes);
sz(sz == 1) = [];
nodes = reshape(nodes, sz);

n1 = nodes(1:end-1,1:end-1);
n2 = nodes(2:end,1:end-1);
n3 = nodes(1:end-1,2:end);
n4 = nodes(2:end,2:end);



f = [n1(:),n2(:),n3(:),n4(:)];

f = f(:,[2 3 1 4]);
end


function f = flipFace(f)
    f = f(:,[1 3 2 4]);   % swap orientation of nodes in faces
end

