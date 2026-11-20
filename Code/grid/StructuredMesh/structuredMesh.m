function grid = structuredMesh(varargin)
%STRUCTUREDMESH Build a Cartesian grid of lines, surfaces, or volumes.
%
%  structuredMesh(x), structuredMesh(x,y), structuredMesh(x,y,z)
%  use prescribed node coordinates. A scalar coordinate fixes that axis.
%
%  structuredMesh(nx,Lx), structuredMesh(nx,ny,Lx,Ly), and
%  structuredMesh(nx,ny,nz,Lx,Ly,Lz) use cell counts and breakpoints.
%  For each varying axis, numel(L) = numel(n)+1. A scalar count of zero
%  fixes an axis at the corresponding scalar coordinate.
%
%  Two positional inputs with a scalar first input use the count form.
%  Use 'input','coordinates' to interpret both inputs as coordinates.
%
%  'triangle',true splits surface quads into triangles.
%  'tetra',true splits volume hexes into tetrahedra.
%
%  grid.nCoord counts coordinate axes and grid.nDim counts varying axes.
%  Line boundary tags 1 and 2 mark the low and high endpoints.
%  Surface tags 1 through 4 mark the low/high second varying axis, then
%  the low/high first varying axis. Volume tags 1 through 6 mark the
%  bottom, top, south, north, west, and east faces.

% Count the positional arguments before the name-value options.
nPos = 0;
for k = 1:numel(varargin)
  if ischar(varargin{k}) || isstring(varargin{k})
    break
  end
  nPos = nPos + 1;
end

options = readInput(struct('tetra',false,'triangle',false,'input',missing), ...
  varargin{nPos+1:end});
inputMode = 'auto';

if ~ismissing(options.input,'missing')
  inputMode = validatestring(options.input,{'auto','coordinates'});
end

countForm = any(nPos == [4,6]) || ...
  (nPos == 2 && isnumeric(varargin{1}) && isscalar(varargin{1}));
axisNames = 'xyz';
if strcmp(inputMode,'auto') && countForm
  nCoord = nPos/2;
  axes = cell(1,nCoord);
  for k = 1:nCoord
    axes{k} = buildStructuredAxis(varargin{k},varargin{nCoord+k},axisNames(k));
  end
elseif ismember(nPos,1:3)
  nCoord = nPos;
  axes = varargin(1:nPos);
else
  error('structuredMesh:InvalidInput', ...
    'Supply one to three coordinate vectors or matching cell counts and breakpoints.');
end

for k = 1:nCoord
  validateattributes(axes{k},{'numeric'},{'vector','nonempty','real','finite'});
  axes{k} = axes{k}(:);
  if any(diff(axes{k}) <= 0)
    error('structuredMesh:InvalidCoordinates', ...
      'Coordinate vectors must be strictly increasing.');
  end
end

axisSizes = cellfun(@numel,axes);
varying = axisSizes(axisSizes > 1);
nDim = numel(varying);
if nDim == 0
  error('structuredMesh:InvalidCoordinates','At least one axis must vary.');
end

tetra = logical(options.tetra);
triangle = logical(options.triangle);
if nDim ~= 3 && tetra
  error('structuredMesh:InvalidInput', ...
    'The tetra option is only available for 3D meshes.');
end
if nDim ~= 2 && triangle
  error('structuredMesh:InvalidInput', ...
    'The triangle option is only available for 2D meshes.');
end

% Coordinates and structured node numbering

components = cell(1,nCoord);
[components{:}] = ndgrid(axes{:});
coord = zeros(prod(axisSizes),nCoord);
for k = 1:nCoord
  coord(:,k) = components{k}(:);
end

if nDim == 1
  nid = (1:varying).';
  topol = [nid(1:end-1),nid(2:end)];
  surf = nid([1,end]);
  surfTag = [1;2];
  cellVTKType = double(VTKType.Line);
  cellNumVerts = 2;
  surfaceVTKType = double(VTKType.Vertex);
  surfaceNumVerts = 1;
elseif nDim == 2
  nid = reshape(1:prod(axisSizes),varying);

  n1 = nid(1:end-1,1:end-1);
  n2 = nid(2:end,  1:end-1);
  n3 = nid(1:end-1,2:end);
  n4 = nid(2:end,  2:end);

  topol = [n1(:),n2(:),n4(:),n3(:)];

  % Boundary lines

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

  nid = reshape(1:prod(axisSizes),varying);
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

  % Select hexahedral or tetrahedral topology

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

% Construct grid

grid = StructuredGrid();
grid.nDim = nDim;
grid.coordinates = coord;
grid.nCoord = size(coord,2);

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

% Mark the grid as structured.
grid.setStructured;

end


function x = buildStructuredAxis(n,L,axisName)

validateattributes(n,{'numeric'},{'vector','nonempty','integer','nonnegative','finite'});
validateattributes(L,{'numeric'},{'vector','nonempty','real','finite'});
n = n(:);
L = L(:);
if isscalar(n) && n == 0 && isscalar(L)
  x = L;
  return
end

if numel(L) ~= numel(n) + 1
  error('structuredMesh:InvalidAxisInput', ...
    'For axis %s, numel(L) must be numel(n)+1.',axisName);
end

if any(n == 0)
  error('structuredMesh:InvalidAxisInput', ...
    'For axis %s, a zero cell count requires a scalar coordinate.',axisName);
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
