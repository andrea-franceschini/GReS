function [matrix, fractureGrids] = createFracturedMesh(bounds, fractures, options)
%CREATEFRACTUREDMESH Mesh a box and split its faces along planar fractures.
%   bounds is [xmin xmax; ymin ymax; zmin zmax]. Each entry of fractures is
%   an ordered N-by-3 perimeter or a struct with center, normal and extents.
%   A single perimeter or a struct array can also be passed directly.
%   extents contains the two full side lengths of a centered rectangle.
%   Its first axis is the least aligned Cartesian axis projected onto the
%   plane (ties use x, then y, then z). Its second axis is normal cross first.
%   matrixMeshSize is the target away from fractures, including the outer
%   boundary. fractureMeshSize is the target near the fractures. If only
%   one is set, it applies everywhere. Omit both for Gmsh defaults.
%   The legacy meshSize option sets one global size.
%   transitionFactor (default 2) multiplies the size difference to set the
%   width of the grading region. Larger values give a more gradual change.
%   fractureGrids{k}.cells.matrixFaces has two matrix.faces indices per
%   triangle. The first face normal follows the perimeter winding or input
%   normal, and the second face normal points in the opposite direction.
%   The mapping is tied to the returned matrix topology.
arguments
  bounds (3,2) double {mustBeReal, mustBeFinite}
  fractures
  options.matrixMeshSize (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 0
  options.fractureMeshSize (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 0
  options.transitionFactor (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 2
  options.meshSize (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 0
end

assert(all(bounds(:,1) < bounds(:,2)), 'GReS:GmshBounds', ...
  'Each box lower bound must be smaller than its upper bound.');
assert(exist('gmshMex','file') == 3, 'GReS:GmshUnavailable', ...
  'Run setupGmsh and compileAll(''gmsh'') to enable fracture meshing.');

if isstruct(fractures)
  fractures = num2cell(fractures);
elseif isnumeric(fractures)
  if isempty(fractures)
    fractures = {};
  else
    fractures = {fractures};
  end
end
assert(iscell(fractures), 'GReS:GmshFractureInput', ...
  'Expected a perimeter, fracture struct, or cell array of fracture definitions.');
polygons = cell(size(fractures));
normals = zeros(numel(fractures),3);
for k = 1:numel(fractures)
  [polygons{k},normals(k,:)] = fracturePerimeter(fractures{k}, k);
end
matrixSize = options.matrixMeshSize;
fractureSize = options.fractureMeshSize;
if options.meshSize > 0
  assert(matrixSize == 0 && fractureSize == 0, 'GReS:GmshMeshSize', ...
    'Use either meshSize or matrixMeshSize and fractureMeshSize.');
  matrixSize = options.meshSize;
end
[coordinates, tetrahedra, triangles] = gmshMex(bounds, polygons, ...
  matrixSize, fractureSize, options.transitionFactor);

matrix = Grid();
matrix.coordinates = coordinates;
matrix.cells.connectivity = tetrahedra;
matrix.cells.tag = ones(size(tetrahedra,1),1);
matrix.cells.numVerts = 4 * matrix.cells.tag;
matrix.cells.VTKType = repmat(double(VTKType.Tetra),size(tetrahedra,1),1);
matrix.processGeometry();

fractureGrids = cell(size(fractures));
if isempty(fractures)
  return
end

matrixFaces = matrix.faces.connectivity.toMatrix();
faceIds = cell(size(triangles));
for k = 1:numel(triangles)
  [matched, faceIds{k}] = ismember(sort(triangles{k},2),matrixFaces,'rows');
  assert(all(matched) && all(~matrix.faces.isBoundary(faceIds{k})), ...
    'GReS:GmshFractureTopology', 'Fracture %d must match internal matrix faces.',k);
end

ids = vertcat(faceIds{:});
assert(numel(unique(ids)) == numel(ids), 'GReS:GmshFractureOverlap', ...
  'Fracture polygons must not overlap over a finite area.');
splitIds = splitMatrixFaces(matrix, ids, matrixFaces);
first = 1;
for k = 1:numel(triangles)
  last = first + size(triangles{k},1) - 1;
  pairs = splitIds(first:last,:);
  swap = matrix.faces.normal(pairs(:,1),:) * normals(k,:).' < 0;
  pairs(swap,:) = pairs(swap,[2 1]);
  fractureGrids{k} = FractureGrid(matrix, triangles{k}, pairs, k);
  first = last + 1;
end
end

function [polygon,normal] = fracturePerimeter(fracture, k)
if isstruct(fracture)
  assert(isscalar(fracture) && all(isfield(fracture,{'center','normal','extents'})), ...
    'GReS:GmshFractureInput', 'Fracture %d needs center, normal and extents.',k);
  validateattributes(fracture.center,{'numeric'},{'real','finite','vector','numel',3});
  validateattributes(fracture.normal,{'numeric'},{'real','finite','vector','numel',3});
  validateattributes(fracture.extents,{'numeric'},{'real','finite','vector','numel',2,'positive'});
  center = double(fracture.center(:)).';
  normal = double(fracture.normal(:)).';
  assert(any(normal), 'GReS:GmshFractureInput', 'Fracture %d has a zero normal.',k);
  normal = normal / max(abs(normal));
  normal = normal / norm(normal);
  [~,axis] = min(abs(normal));
  u = zeros(1,3);
  u(axis) = 1;
  u = u - dot(u,normal)*normal;
  u = u / norm(u);
  v = cross(normal,u);
  extents = double(fracture.extents);
  polygon = center + [-1;1;1;-1] * (extents(1)/2*u) + ...
                     [-1;-1;1;1] * (extents(2)/2*v);
else
  validateattributes(fracture,{'numeric'},{'real','finite','2d','ncols',3});
  polygon = full(double(fracture));
end
if size(polygon,1) > 1 && isequal(polygon(1,:),polygon(end,:))
  polygon(end,:) = [];
end
assert(size(polygon,1) >= 3 && size(unique(polygon,'rows'),1) == size(polygon,1), ...
  'GReS:GmshFractureInput', 'Fracture %d needs at least three distinct vertices.',k);
centered = polygon - mean(polygon,1);
span = svd(centered,0);
tol = 1e-10 * span(1);
assert(span(2) > tol && span(3) <= tol, 'GReS:GmshFractureInput', ...
  'Fracture %d must be planar and nondegenerate.',k);
if ~isstruct(fracture)
  normal = sum(cross(centered,circshift(centered,-1,1),2),1);
  assert(norm(normal) > tol * span(2), 'GReS:GmshFractureInput', ...
    'Fracture %d must enclose a nonzero area.',k);
  normal = normal / norm(normal);
end
end

function pairs = splitMatrixFaces(matrix, ids, connectivity)
faces = matrix.faces;
otherCells = faces.neighbors(ids,2);
newIds = faces.num + (1:numel(ids)).';
pairs = [ids,newIds];

% Redirect only the second cell to the duplicate face. Keep local face IDs.
cellFaces = matrix.cells.cells2faces.toMatrix();
[rows,localFaces] = find(cellFaces(otherCells,:) == ids);
assert(numel(rows) == numel(ids), 'GReS:GmshFractureTopology', ...
  'Could not locate both matrix sides of each fracture.');
cellFaces(sub2ind(size(cellFaces),otherCells(rows),localFaces)) = newIds(rows);
matrix.cells.cells2faces = ArrayOfArrays(cellFaces);

faces.num = faces.num + numel(ids);
faces.connectivity = ArrayOfArrays([connectivity; connectivity(ids,:)]);
faces.neighbors(ids,2) = 0;
faces.neighbors(newIds,:) = [otherCells,zeros(numel(ids),1)];
faces.isBoundary(pairs(:)) = true;
faces.normal(newIds,:) = -faces.normal(ids,:);
faces.center(newIds,:) = faces.center(ids,:);
faces.area(newIds,1) = faces.area(ids);
faces.numVerts(newIds,1) = faces.numVerts(ids);
faces.VTKType(newIds,1) = faces.VTKType(ids);
matrix.faces = faces;
end
