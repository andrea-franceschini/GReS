function processGeometry(grid,preserveParentGeometry)
% Populate cell measures and codimension-one geometry for 2D and 3D grids.
% An extracted boundary patch inherits its parent face geometry.
fromParent = nargin > 1 && preserveParentGeometry;
if fromParent
  parentGeometry = grid.cells;
end
initializeGrid(grid);
assert(grid.cells.num > 0 && ismember(grid.nDim,[2,3]), ...
  'GReS:GridDimension', 'Expected a populated 2D or 3D grid.');

if fromParent
  grid.cells.normal = parentGeometry.normal;
elseif grid.nDim == 2
  previousNormals = [];
  if isfield(grid.cells,'normal')
    previousNormals = grid.cells.normal;
  end
  grid.cells.normal = zeros(grid.cells.num,3);
  for vtkId = grid.cells.vtkTypes
    ids = grid.getCellsByVTKId(vtkId);
    nodes = grid.getCellNodes(ids);
    nodes = nodes(:,1:grid.cells.numVerts(ids(1)));
    if vtkId == VTKType.Quad9
      nodes = nodes(:,1:4);
    end
    xyz = grid.coordinates(reshape(nodes.',[],1),:);
    if grid.nCoord == 2
      xyz(:,3) = 0;
    end
    [~,~,normal] = computePolygonGeometry(xyz,repmat(size(nodes,2),numel(ids),1));
    if ~isempty(previousNormals)
      flip = dot(normal,previousNormals(ids,:),2) < 0;
      normal(flip,:) = -normal(flip,:);
    end
    grid.cells.normal(ids,:) = normal;
  end
end

grid.processFaces();
if fromParent
  grid.cells.volume = parentGeometry.area;
  grid.cells.center = parentGeometry.center;
else
  for vtkId = grid.cells.vtkTypes
    ids = grid.getCellsByVTKId(vtkId);
    elem = FiniteElementType.create(vtkId,grid);
    [grid.cells.volume(ids),grid.cells.center(ids,:)] = elem.getSizeAndCentroid(ids);
  end
end
if grid.nDim == 2
  grid.cells.area = grid.cells.volume;
end

neighbors = grid.faces.neighbors(:,1);
outward = grid.faces.center - grid.cells.center(neighbors,:);
flip = dot(outward,grid.faces.normal,2) < 0;
grid.faces.normal(flip,:) = -grid.faces.normal(flip,:);

faceIds = grid.surfaces.faceId;
grid.surfaces.area = grid.faces.area(faceIds);
grid.surfaces.center = grid.faces.center(faceIds,:);
grid.surfaces.normal = grid.faces.normal(faceIds,:);
grid.isProcessed = true;
end
