function grid = structuredMesh(nx,ny,nz,Lx,Ly,Lz)
% STRUCTUREDMESH  Build a Cartesian structured hexahedral mesh.
%
%   GRID = STRUCTUREDMESH(nx, ny, nz, Lx, Ly, Lz) creates a structured
%   hexahedral mesh over the axis-aligned box defined by the breakpoint
%   arrays Lx, Ly, Lz.
%
%   Here:
%     - nx, ny, nz are integer arrays
%     - Lx, Ly, Lz are coordinate arrays with one element more
%     - nx(i) is the number of cells between Lx(i) and Lx(i+1)
%       and similarly for y and z
%
%   Example:
%     % x: 4 cells on [0,1], 8 cells on [1,3]
%     % y: 3 cells on [0,2]
%     % z: 2 cells on [0,0.5], 1 cell on [0.5,1]
%     grid = structuredMesh([4 8], [3], [2 1], ...
%                           [0 1 3], [0 2], [0 0.5 1]);
%
%   GRID = STRUCTUREDMESH(x, y, z) creates a mesh whose node positions are
%   given explicitly by the coordinate vectors x, y, z.
%
% -------------------------------------------------------------------------
% INPUTS (piecewise-uniform spacing)
%   nx, ny, nz  - number of cells in each interval      [integer arrays]
%   Lx, Ly, Lz  - interval breakpoints                  [double arrays]
%                 with numel(L*) = numel(n*) + 1
%
% INPUTS (prescribed node coordinates)
%   x, y, z     - node coordinate vectors               [double arrays]
%
% -------------------------------------------------------------------------
% OUTPUT
%   grid  - grid object with processed geometry

if nargin == 6
  x = buildStructuredAxis(nx, Lx, 'x');
  y = buildStructuredAxis(ny, Ly, 'y');
  z = buildStructuredAxis(nz, Lz, 'z');

  nnx = numel(x);
  nny = numel(y);
  nnz = numel(z);

elseif nargin == 3
  x = nx(:);
  y = ny(:);
  z = nz(:);

  nnx = numel(x);
  nny = numel(y);
  nnz = numel(z);

else
  error('structuredMesh:InvalidInput', ...
      'Use either structuredMesh(nx,ny,nz,Lx,Ly,Lz) or structuredMesh(x,y,z).');
end

[X,Y,Z] = ndgrid(x,y,z);
coord = [X(:), Y(:), Z(:)];

% linear node numbering
nid = reshape(1:nnx*nny*nnz, nnx, nny, nnz);

n1 = nid(1:end-1,1:end-1,1:end-1);
n2 = nid(2:end,  1:end-1,1:end-1);
n3 = nid(1:end-1,2:end,  1:end-1);
n4 = nid(2:end,  2:end,  1:end-1);
n5 = nid(1:end-1,1:end-1,2:end);
n6 = nid(2:end,  1:end-1,2:end);
n7 = nid(1:end-1,2:end,  2:end);
n8 = nid(2:end,  2:end,  2:end);

topol = [n1(:),n2(:),n4(:),n3(:),n5(:),n6(:),n8(:),n7(:)];

% boundary surfaces

% bottom: tag 1
n = nid(:,:,1);
fb = flipFace(FaceGrid(n));

% top: tag 2
n = nid(:,:,end);
ft = FaceGrid(n);

% south: tag 3
n = nid(:,1,:);
fs = FaceGrid(n);

% north: tag 4
n = nid(:,end,:);
fn = flipFace(FaceGrid(n));

% west: tag 5
n = nid(1,:,:);
fw = flipFace(FaceGrid(n));

% east: tag 6
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
grid.cells.VTKType      = 12 * ones(nC,1);
grid.cells.tag          = ones(nC,1);
grid.cells.numVerts     = 8  * ones(nC,1);

nS = size(surf,1);
grid.surfaces.connectivity = surf;
grid.surfaces.VTKType      = 9 * ones(nS,1);
grid.surfaces.tag          = surfTag;
grid.surfaces.numVerts     = 4 * ones(nS,1);

% finalize geometry
grid.setStructured;
processGeometry(grid);

end


function x = buildStructuredAxis(n, L, axisName)

n = n(:);
L = L(:);

if numel(L) ~= numel(n) + 1
  error('structuredMesh:InvalidAxisInput', ...
      'For axis %s, numel(L) must be numel(n)+1.', axisName);
end

if any(n <= 0) || any(abs(n - round(n)) > 0)
  error('structuredMesh:InvalidAxisInput', ...
      'For axis %s, cell counts must be positive integers.', axisName);
end

if any(diff(L) <= 0)
  error('structuredMesh:InvalidAxisInput', ...
      'For axis %s, breakpoints must be strictly increasing.', axisName);
end

x = zeros(sum(n) + 1, 1);
k = 1;

for i = 1:numel(n)
  xi = linspace(L(i), L(i+1), n(i)+1).';
  if i < numel(n)
    xi = xi(1:end-1);
  end
  nk = numel(xi);
  x(k:k+nk-1) = xi;
  k = k + nk;
end

if k <= numel(x)
  x(k:end) = L(end);
end

end


function f = FaceGrid(nodes)

sz = size(nodes);
sz(sz == 1) = [];
nodes = reshape(nodes, sz);

n1 = nodes(1:end-1,1:end-1);
n2 = nodes(2:end,  1:end-1);
n3 = nodes(1:end-1,2:end);
n4 = nodes(2:end,  2:end);

f = [n1(:), n2(:), n4(:), n3(:)];

end


function f = flipFace(f)
f = f(:, [1 4 3 2]);
end