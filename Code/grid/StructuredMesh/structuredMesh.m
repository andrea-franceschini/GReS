function grid = structuredMesh(nx,ny,varargin)
% STRUCTUREDMESH Build a 2D or 3D Cartesian structured mesh.
%
%   GRID = STRUCTUREDMESH(nx,ny,Lx,Ly)
%   creates a structured 2D quadrilateral mesh.
%
%   GRID = STRUCTUREDMESH(x,y)
%   creates a quadrilateral mesh using prescribed node coordinates.
%
%   GRID = STRUCTUREDMESH(nx,ny,nz,Lx,Ly,Lz)
%   creates a structured 3D hexahedral mesh.
%
%   GRID = STRUCTUREDMESH(x,y,z)
%   creates a hexahedral mesh using prescribed node coordinates.
%
%   GRID = STRUCTUREDMESH(...,'triangle',1)
%   divides each 2D quadrilateral into two triangles.
%
%   GRID = STRUCTUREDMESH(...,'tetra',1)
%   divides each 3D hexahedron into six tetrahedra.
%
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
%   x, y, z     - strictly increasing node coordinate vectors
%
% OPTIONAL INPUT
%
%   'triangle'  - false/0: quadrilateral mesh, default
%                 true/1: triangular mesh (2D only)
%
%   'tetra'     - false/0: hexahedral mesh, default
%                 true/1: tetrahedral mesh (3D only)
%
% OUTPUT
%
%   grid        - processed grid, with grid.nDim = 2 or 3
%
%   2D coordinates are stored as [x,y,0]. Cells are quads or triangles;
%   boundary lines are stored in grid.surfaces (VTK_LINE = 3).
%   2D tags: bottom (y min) = 1, top (y max) = 2,
%            left (x min) = 3, right (x max) = 4.
%   3D tags remain: bottom = 1, top = 2, south = 3, north = 4,
%                    west = 5, east = 6.

% Count the positional arguments before the name-value options.
nPos = 0;
for k = 1:numel(varargin)
  if ischar(varargin{k}) || isstring(varargin{k})
    break
  end
  nPos = nPos + 1;
end

switch nPos
  case 0
    % structuredMesh(x,y,...)
    nDim = 2;
    x = nx(:);
    y = ny(:);
    z = 0;

  case 1
    % structuredMesh(x,y,z,...)
    nDim = 3;
    x = nx(:);
    y = ny(:);
    z = varargin{1}(:);

  case 2
    % structuredMesh(nx,ny,Lx,Ly,...)
    nDim = 2;
    x = buildStructuredAxis(nx,varargin{1},'x');
    y = buildStructuredAxis(ny,varargin{2},'y');
    z = 0;

  case 4
    % structuredMesh(nx,ny,nz,Lx,Ly,Lz,...)
    nDim = 3;
    nz = varargin{1};
    x = buildStructuredAxis(nx,varargin{2},'x');
    y = buildStructuredAxis(ny,varargin{3},'y');
    z = buildStructuredAxis(nz,varargin{4},'z');

  otherwise
    error('structuredMesh:InvalidInput', ...
      ['Use structuredMesh(x,y), structuredMesh(x,y,z), ', ...
       'structuredMesh(nx,ny,Lx,Ly), or ', ...
       'structuredMesh(nx,ny,nz,Lx,Ly,Lz).']);
end

inputOptions = varargin(nPos+1:end);

% Parse optional arguments.
options = struct('tetra',false,'triangle',false);
options = readInput(options,inputOptions{:});
tetra = logical(options.tetra);
triangle = logical(options.triangle);

if nDim == 2 && tetra
  error('structuredMesh:InvalidInput', ...
    'The tetra option is only available for 3D meshes.');
end

if nDim == 3 && triangle
  error('structuredMesh:InvalidInput', ...
    'The triangle option is only available for 2D meshes.');
end

nnx = numel(x);
nny = numel(y);
nnz = numel(z);

if nnx < 2 || nny < 2 || (nDim == 3 && nnz < 2)
  error('structuredMesh:InvalidCoordinates', ...
    'Each active coordinate vector must contain at least two entries.');
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

if nDim == 2

  n1 = nid(1:end-1,1:end-1);
  n2 = nid(2:end,  1:end-1);
  n3 = nid(1:end-1,2:end);
  n4 = nid(2:end,  2:end);

  topol = [n1(:),n2(:),n4(:),n3(:)];

  %% Boundary lines

  % Bottom (constant y): tag 1
  n = nid(:,1);
  fs = [n(1:end-1),n(2:end)];

  % Top (constant y): tag 2
  n = nid(:,end);
  fn = [n(2:end),n(1:end-1)];

  % Left (constant x): tag 3
  n = nid(1,:).';
  fw = [n(2:end),n(1:end-1)];

  % Right (constant x): tag 4
  n = nid(end,:).';
  fe = [n(1:end-1),n(2:end)];

  surf = [fs;fn;fw;fe];

  surfTag = [
    ones(size(fs,1),1)
    2*ones(size(fn,1),1)
    3*ones(size(fw,1),1)
    4*ones(size(fe,1),1)
  ];

  if triangle

    topol = QuadToTriangle(topol);
    cellVTKType = 5; % VTK_TRIANGLE
    cellNumVerts = 3;

  else

    cellVTKType = 9; % VTK_QUAD
    cellNumVerts = 4;

  end

  surfaceVTKType = 3; % VTK_LINE
  surfaceNumVerts = 2;

else

  n1 = nid(1:end-1,1:end-1,1:end-1);
  n2 = nid(2:end,  1:end-1,1:end-1);
  n3 = nid(1:end-1,2:end,  1:end-1);
  n4 = nid(2:end,  2:end,  1:end-1);
  n5 = nid(1:end-1,1:end-1,2:end);
  n6 = nid(2:end,  1:end-1,2:end);
  n7 = nid(1:end-1,2:end,  2:end);
  n8 = nid(2:end,  2:end,  2:end);

  % Hexahedral local ordering:
  % local [1 2 3 4 5 6 7 8] = [n1 n2 n4 n3 n5 n6 n8 n7]
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

    topol = HexaToTetra(topolHex);
    surf = QuadToTriangle(surfQuad);
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

end

%% Construct grid

grid = Grid();
grid.nDim = nDim;
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
% All tetrahedra share the body diagonal between local vertices 1 and 7.

nHexa = size(hexa,1);

tetra = zeros(6*nHexa,4,'like',hexa);

tetra(1:6:end,:) = hexa(:,[1 2 3 7]);
tetra(2:6:end,:) = hexa(:,[1 3 4 7]);
tetra(3:6:end,:) = hexa(:,[1 4 8 7]);
tetra(4:6:end,:) = hexa(:,[1 8 5 7]);
tetra(5:6:end,:) = hexa(:,[1 5 6 7]);
tetra(6:6:end,:) = hexa(:,[1 6 2 7]);

end


function triangle = QuadToTriangle(quad)
% QUADTOTRIANGLE Divide every oriented quadrilateral into two triangles.

nQuad = size(quad,1);

triangle = zeros(2*nQuad,3,'like',quad);

triangle(1:2:end,:) = quad(:,[1 2 3]);
triangle(2:2:end,:) = quad(:,[1 3 4]);

end
