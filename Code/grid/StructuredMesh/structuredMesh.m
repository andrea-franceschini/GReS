function grid = structuredMesh(nx,ny,nz,varargin)
% STRUCTUREDMESH Build a Cartesian structured hexahedral or tetrahedral mesh.
%
%   GRID = STRUCTUREDMESH(nx,ny,nz,Lx,Ly,Lz)
%   creates a structured hexahedral mesh.
%
%   GRID = STRUCTUREDMESH(nx,ny,nz,Lx,Ly,Lz,'tetra',1)
%   creates the same structured grid, with every hexahedron divided into
%   six tetrahedra.
%
%   GRID = STRUCTUREDMESH(x,y,z)
%   creates a hexahedral mesh using explicitly prescribed node coordinates.
%
%   GRID = STRUCTUREDMESH(x,y,z,'tetra',1)
%   creates the corresponding tetrahedral mesh.
%
% -------------------------------------------------------------------------
% INPUTS, piecewise-uniform spacing
%
%   nx, ny, nz  - number of cells in each interval
%   Lx, Ly, Lz  - interval breakpoints, with:
%
%                   numel(Lx) = numel(nx) + 1
%
%                 and similarly for y and z.
%
% INPUTS, prescribed coordinates
%
%   x, y, z     - node coordinate vectors
%
% OPTIONAL INPUT
%
%   'tetra'     - false/0: hexahedral mesh, default
%                 true/1:  tetrahedral mesh
%
% -------------------------------------------------------------------------
% OUTPUT
%
%   grid        - processed three-dimensional grid

% Determine whether the first three optional inputs are the breakpoint
% arrays or whether varargin contains only name-value arguments.
if isempty(varargin) || ischar(varargin{1}) || isstring(varargin{1})

  % Explicit coordinate-vector input:
  %
  % structuredMesh(x,y,z,...)

  x = nx(:);
  y = ny(:);
  z = nz(:);

  inputOptions = varargin;

else

  % Piecewise-uniform input:
  %
  % structuredMesh(nx,ny,nz,Lx,Ly,Lz,...)

  if numel(varargin) < 3
    error('structuredMesh:InvalidInput', ...
      ['Use structuredMesh(nx,ny,nz,Lx,Ly,Lz) or ', ...
       'structuredMesh(x,y,z).']);
  end

  Lx = varargin{1};
  Ly = varargin{2};
  Lz = varargin{3};

  x = buildStructuredAxis(nx,Lx,'x');
  y = buildStructuredAxis(ny,Ly,'y');
  z = buildStructuredAxis(nz,Lz,'z');

  inputOptions = varargin(4:end);

end

% Parse optional arguments.
options = struct('tetra',0);
options = readInput(options,inputOptions{:});
tetra = logical(options.tetra);

nnx = numel(x);
nny = numel(y);
nnz = numel(z);

if nnx < 2 || nny < 2 || nnz < 2
  error('structuredMesh:InvalidCoordinates', ...
    'Each coordinate vector must contain at least two entries.');
end

if any(diff(x) <= 0) || any(diff(y) <= 0) || any(diff(z) <= 0)
  error('structuredMesh:InvalidCoordinates', ...
    'Coordinate vectors must be strictly increasing.');
end

%% Coordinates

[X,Y,Z] = ndgrid(x,y,z);
coord = [X(:),Y(:),Z(:)];

%% Structured node numbering

nid = reshape(1:nnx*nny*nnz,nnx,nny,nnz);

n1 = nid(1:end-1,1:end-1,1:end-1);
n2 = nid(2:end,  1:end-1,1:end-1);
n3 = nid(1:end-1,2:end,  1:end-1);
n4 = nid(2:end,  2:end,  1:end-1);
n5 = nid(1:end-1,1:end-1,2:end);
n6 = nid(2:end,  1:end-1,2:end);
n7 = nid(1:end-1,2:end,  2:end);
n8 = nid(2:end,  2:end,  2:end);

% Hexahedral local ordering:
%
%       8---------7
%      /|        /|
%     5---------6 |
%     | |       | |
%     | 4-------|-3
%     |/        |/
%     1---------2
%
% In terms of the structured arrays:
%
%   local 1 = n1
%   local 2 = n2
%   local 3 = n4
%   local 4 = n3
%   local 5 = n5
%   local 6 = n6
%   local 7 = n8
%   local 8 = n7

topolHex = [
  n1(:), ...
  n2(:), ...
  n4(:), ...
  n3(:), ...
  n5(:), ...
  n6(:), ...
  n8(:), ...
  n7(:)
];

%% Boundary quadrilateral surfaces

% Bottom: tag 1
n = nid(:,:,1);
fb = flipFace(FaceGrid(n));

% Top: tag 2
n = nid(:,:,end);
ft = FaceGrid(n);

% South: tag 3
n = nid(:,1,:);
fs = FaceGrid(n);

% North: tag 4
n = nid(:,end,:);
fn = flipFace(FaceGrid(n));

% West: tag 5
n = nid(1,:,:);
fw = flipFace(FaceGrid(n));

% East: tag 6
n = nid(end,:,:);
fe = FaceGrid(n);

surfQuad = [fb;ft;fs;fn;fw;fe];

surfTagQuad = [
  ones(size(fb,1),1)
  2*ones(size(ft,1),1)
  3*ones(size(fs,1),1)
  4*ones(size(fn,1),1)
  5*ones(size(fw,1),1)
  6*ones(size(fe,1),1)
];

%% Select hexahedral or tetrahedral topology

if tetra

  % Split every hexahedron into six tetrahedra around the common body
  % diagonal connecting local hexahedral vertices 1 and 7.
  topol = HexaToTetra(topolHex);

  % Split each oriented quadrilateral boundary face into two triangles.
  surf = QuadToTriangle(surfQuad);

  % QuadToTriangle places the two triangles associated with each
  % quadrilateral consecutively.
  surfTag = repelem(surfTagQuad,2,1);

  cellVTKType = 10; % VTK_TETRA
  cellNumVerts = 4;

  surfaceVTKType = 5; % VTK_TRIANGLE
  surfaceNumVerts = 3;

else

  topol = topolHex;
  surf = surfQuad;
  surfTag = surfTagQuad;

  cellVTKType = 12; % VTK_HEXAHEDRON
  cellNumVerts = 8;

  surfaceVTKType = 9; % VTK_QUAD
  surfaceNumVerts = 4;

end

%% Construct grid

grid = Grid();
grid.nDim = 3;
grid.coordinates = coord;

nC = size(topol,1);

grid.cells.connectivity = topol;
grid.cells.VTKType = cellVTKType*ones(nC,1);
grid.cells.tag = ones(nC,1);
grid.cells.numVerts = cellNumVerts*ones(nC,1);

nS = size(surf,1);

grid.surfaces.connectivity = surf;
grid.surfaces.VTKType = surfaceVTKType*ones(nS,1);
grid.surfaces.tag = surfTag;
grid.surfaces.numVerts = surfaceNumVerts*ones(nS,1);

% Finalize geometry.
grid.setStructured;

end


function x = buildStructuredAxis(n,L,axisName)

n = n(:);
L = L(:);

if numel(L) ~= numel(n) + 1
  error('structuredMesh:InvalidAxisInput', ...
    'For axis %s, numel(L) must be numel(n)+1.',axisName);
end

if any(n <= 0) || any(abs(n-round(n)) > 0)
  error('structuredMesh:InvalidAxisInput', ...
    'For axis %s, cell counts must be positive integers.',axisName);
end

if any(diff(L) <= 0)
  error('structuredMesh:InvalidAxisInput', ...
    'For axis %s, breakpoints must be strictly increasing.',axisName);
end

x = zeros(sum(n)+1,1);
k = 1;

for i = 1:numel(n)

  xi = linspace(L(i),L(i+1),n(i)+1).';

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
nodes = reshape(nodes,sz);

n1 = nodes(1:end-1,1:end-1);
n2 = nodes(2:end,  1:end-1);
n3 = nodes(1:end-1,2:end);
n4 = nodes(2:end,  2:end);

f = [n1(:),n2(:),n4(:),n3(:)];

end


function f = flipFace(f)

f = f(:,[1 4 3 2]);

end


function tetra = HexaToTetra(hexa)
% HEXATOTETRA Divide each hexahedron into six tetrahedra.
%
% Input:
%
%   hexa  - nHexa x 8
%
% Output:
%
%   tetra - 6*nHexa x 4
%
% All tetrahedra share the hexahedral body diagonal between local
% vertices 1 and 7. This gives a conforming subdivision across adjacent
% structured cells.

nHexa = size(hexa,1);

tetra = zeros(6*nHexa,4,'like',hexa);

% Local tetrahedral connectivity:
%
%   [1 2 3 7]
%   [1 3 4 7]
%   [1 4 8 7]
%   [1 8 5 7]
%   [1 5 6 7]
%   [1 6 2 7]

tetra(1:6:end,:) = hexa(:,[1 2 3 7]);
tetra(2:6:end,:) = hexa(:,[1 3 4 7]);
tetra(3:6:end,:) = hexa(:,[1 4 8 7]);
tetra(4:6:end,:) = hexa(:,[1 8 5 7]);
tetra(5:6:end,:) = hexa(:,[1 5 6 7]);
tetra(6:6:end,:) = hexa(:,[1 6 2 7]);

end


function triangle = QuadToTriangle(quad)
% QUADTOTRIANGLE Divide every oriented quadrilateral into two triangles.
%
% Input:
%
%   quad      - nQuad x 4
%
% Output:
%
%   triangle  - 2*nQuad x 3
%
% The orientation inherited from the quadrilateral is preserved.

nQuad = size(quad,1);

triangle = zeros(2*nQuad,3,'like',quad);

triangle(1:2:end,:) = quad(:,[1 2 3]);
triangle(2:2:end,:) = quad(:,[1 3 4]);

end