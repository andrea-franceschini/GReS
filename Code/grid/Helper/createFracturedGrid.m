function [matrix, fractureGrids] = createFracturedGrid(bounds, fractures, varargin)
%CREATEFRACTUREDGRID Mesh a box and split its faces along planar fractures.
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
%   fracture cell. The first face normal follows the perimeter winding or input
%   normal, and the second face normal points in the opposite direction.
%   structured=true builds a Cartesian grid with coordinate-aligned
%   fracture planes and perimeters, using the smaller supplied mesh size
%   throughout. transitionFactor applies only to unstructured grids.
%   Matrix nodes are shared across intact faces, including internal fracture
%   fronts. Separate rock sectors at fracture interiors and crossings use
%   coincident node copies.


% deal input parameters
default = struct('matrixMeshSize',0,...
                 'fractureMeshSize',0,...
                 'transitionFactor',2.0,...
                 'meshSize',0,...
                 'structured',false);

options = readInput(default,varargin{:});

assert(all(size(bounds) == [3,2]));
assert(all(bounds(:,1) < bounds(:,2)), 'GReS:GmshBounds', ...
  'Each box lower bound must be smaller than its upper bound.');

% collect fractures in a cell array
if isstruct(fractures)
  fractures = num2cell(fractures);
elseif isnumeric(fractures)
  if isempty(fractures)
    fractures = {};
  else
    fractures = {fractures};
  end
end

polygons = cell(size(fractures));
normals = zeros(numel(fractures),3);

% get unified fracture geometry description
% resolve struct fracture as perimeters
for k = 1:numel(fractures)
  [polygons{k},normals(k,:)] = fracturePerimeter(fractures{k}, k);
end

% requested mesh size for matrix and fracture
matrixSize = options.matrixMeshSize;
fractureSize = options.fractureMeshSize;
if options.meshSize > 0
  assert(matrixSize == 0 && fractureSize == 0, 'GReS:GmshMeshSize', ...
    'Use either meshSize or matrixMeshSize and fractureMeshSize.');
  matrixSize = double(options.meshSize);
end

if options.structured
  % attempt to create a structured grid with hexahedra 
  % requires the fracture to be parallel to the coordinate axis
  [matrix,faceIds] = structuredFractureMesh(bounds,polygons,normals, ...
    matrixSize,fractureSize);
else

  % Build the conforming tetrahedral mesh through the individual Gmsh calls.
  [coordinates, matrixCells, fractureCells] = gresgmsh.meshFracturedBox( ...
    bounds, polygons, matrixSize, fractureSize, options.transitionFactor);

  matrix = Grid();
  matrix.coordinates = coordinates;
  matrix.cells.connectivity = matrixCells;
  matrix.cells.tag = ones(size(matrixCells,1),1);
  matrix.cells.numVerts = repmat(size(matrixCells,2),size(matrixCells,1),1);
  matrix.cells.VTKType = repmat(double(VTKType.Tetra),size(matrixCells,1),1);
  matrix.processGeometry();

  % check if fractures and matrix are conforming
  faceIds = cell(size(fractureCells));
  matrixFaces = sort(matrix.faces.connectivity.toMatrix(),2);
  for k = 1:numel(fractureCells)
    width = size(fractureCells{k},2);

    [matched, faceIds{k}] = ismember(sort(fractureCells{k},2), ...
      matrixFaces(:,end-width+1:end),'rows');
    assert(all(matched), 'GReS:GmshFractureTopology', ...
      'Fracture %d must match matrix faces.',k);
  end
end

fractureGrids = cell(size(fractures));
if isempty(fractures)
  return
end

for k = 1:numel(faceIds)
  assert(~isempty(faceIds{k}) && all(~matrix.faces.isBoundary(faceIds{k})), ...
    'GReS:GmshFractureTopology', 'Fracture %d must match internal matrix faces.',k);
end

ids = vertcat(faceIds{:});
assert(numel(unique(ids)) == numel(ids), 'GReS:GmshFractureOverlap', ...
  'Fracture polygons must not overlap over a finite area.');
matrixFaces = matrix.faces.connectivity.toMatrix();

splitIds = splitMatrixFaces(matrix, ids, matrixFaces);
first = 1;
for k = 1:numel(faceIds)
  last = first + numel(faceIds{k}) - 1;
  pairs = splitIds(first:last,:);
  swap = matrix.faces.normal(pairs(:,1),:) * normals(k,:).' < 0;
  pairs(swap,:) = pairs(swap,[2 1]);
  fractureGrids{k} = FractureGrid(matrix, matrixFaces(faceIds{k},:), pairs, k);
  if options.structured
    fractureGrids{k}.setStructured();
  end
  first = last + 1;
end
end

function [matrix,faceIds] = structuredFractureMesh(bounds,polygons,normals,matrixSize,fractureSize)
axes = num2cell(bounds.',1);
normalAxes = zeros(numel(polygons),1);
for k = 1:numel(polygons)
  polygon = polygons{k};
  [~,axis] = max(abs(normals(k,:)));
  tangential = setdiff(1:3,axis);
  edges = circshift(polygon,-1,1) - polygon;
  assert(all(abs(normals(k,tangential)) <= 1e-10) && ...
    all(polygon(:,axis) == polygon(1,axis)) && ...
    all(sum(edges(:,tangential) ~= 0,2) == 1), ...
    'GReS:StructuredFractureAlignment', ...
    'Fracture %d must have a coordinate-aligned plane and perimeter.',k);
  assert(all(all(polygon >= bounds(:,1).' & polygon <= bounds(:,2).')), ...
    'GReS:StructuredFractureBounds','Fracture %d must lie inside the box.',k);
  normalAxes(k) = axis;
  for d = 1:3
    axes{d} = [axes{d};polygon(:,d)];
  end
end

sizes = [matrixSize,fractureSize];
step = min(sizes(sizes > 0));
for d = 1:3
  breaks = unique(axes{d});
  if isempty(step)
    axes{d} = breaks;
    continue
  end
  intervals = cell(numel(breaks)-1,1);
  for j = 1:numel(intervals)
    n = ceil((breaks(j+1)-breaks(j))/step);
    points = linspace(breaks(j),breaks(j+1),n+1).';
    intervals{j} = points(1:end-1);
  end
  axes{d} = [vertcat(intervals{:});breaks(end)];
end
matrix = structuredMesh(axes{:});
matrix.processGeometry();

centers = matrix.faces.center;
faceNodes = matrix.faces.connectivity.toMatrix();
faceIds = cell(size(polygons));
for k = 1:numel(polygons)
  axis = normalAxes(k);
  tangential = setdiff(1:3,axis);
  polygon = polygons{k};
  candidates = find(matrix.coordinates(faceNodes(:,1),axis) == polygon(1,axis) & ...
    abs(matrix.faces.normal(:,axis)) > 0.5);
  inside = inpolygon(centers(candidates,tangential(1)),centers(candidates,tangential(2)), ...
    polygon(:,tangential(1)),polygon(:,tangential(2)));
  faceIds{k} = candidates(inside);
end
end

function [polygon,normal] = fracturePerimeter(fracture, k)

if isstruct(fracture)
  assert(all(isfield(fracture,{'center','normal','extents'})), ...
    'GReS:GmshFractureInput', 'Fracture %d needs center, normal and extents.',k);
  
  center = fracture.center(:).';
  normal = fracture.normal(:).';
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
oldCells = matrix.cells.connectivity;
newCells = splitMatrixNodes(matrix, ids, connectivity);
matrix.cells.connectivity = newCells;
matrix.nNodes = size(matrix.coordinates,1);
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
owners = [faces.neighbors(:,1);otherCells];
faces.connectivity = ArrayOfArrays(remapNodes([connectivity;connectivity(ids,:)], ...
  oldCells(owners,:),newCells(owners,:)));
faces.neighbors(ids,2) = 0;
faces.neighbors(newIds,:) = [otherCells,zeros(numel(ids),1)];
faces.isBoundary(pairs(:)) = true;
faces.normal(newIds,:) = -faces.normal(ids,:);
faces.center(newIds,:) = faces.center(ids,:);
faces.area(newIds,1) = faces.area(ids);
faces.numVerts(newIds,1) = faces.numVerts(ids);
faces.VTKType(newIds,1) = faces.VTKType(ids);
matrix.faces = faces;
owners = faces.neighbors(matrix.surfaces.faceId,1);
surfaceNodes = matrix.getSurfNodes((1:matrix.surfaces.num).');
matrix.surfaces.connectivity = remapNodes(surfaceNodes,oldCells(owners,:),newCells(owners,:));
end

function cells = splitMatrixNodes(matrix, ids, connectivity)
% Give each connected component of the node's incident cells one node copy.
% Connections cross intact faces only. Cells remain connected around an

% internal fracture front, so its nodes stay shared.
cells = matrix.cells.connectivity;
neighbors = matrix.faces.neighbors;
isSplit = false(matrix.faces.num,1);
isSplit(ids) = true;
splitNodes = unique(connectivity(ids,:));
splitNodes(splitNodes == 0) = [];
[isSplitNode,nodeIndex] = ismember(connectivity,splitNodes);
[faceRows,~] = find(isSplitNode);
incidentFaces = accumarray(nodeIndex(isSplitNode),faceRows,[numel(splitNodes),1], ...
  @(rows) {rows});
coordinates = matrix.coordinates;

for k = 1:numel(splitNodes)
  node = splitNodes(k);
  faceRows = incidentFaces{k};
  owners = neighbors(faceRows,:);
  [incidentCells,~,local] = unique(owners(:));
  local = reshape(local,size(owners));
  internal = owners(:,2) ~= 0;
  ties = local(internal & ~isSplit(faceRows),:);
  groups = (1:numel(incidentCells)).';

  for j = 1:size(ties,1)
    a = groups(ties(j,1));
    b = groups(ties(j,2));
    groups(groups == b) = a;
  end

  [~,~,group] = unique(groups(incidentCells > 0),'stable');
  newIds = size(coordinates,1) + (1:max(group)-1);
  coordinates(newIds,:) = repmat(coordinates(node,:),numel(newIds),1);
  rows = incidentCells(incidentCells > 0);
  copies = [node;newIds(:)];
  cells(rows,:) = remapNodes(cells(rows,:),repmat(node,numel(rows),1), ...
    copies(group));
end

matrix.coordinates = coordinates;
end

function nodes = remapNodes(nodes,oldCells,newCells)
original = nodes;
  for j = 1:size(oldCells,2)
    [rows,cols] = find(original == oldCells(:,j));
    nodes(sub2ind(size(nodes),rows,cols)) = newCells(rows,j);
  end
end
