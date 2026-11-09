function processFaces(grid)
% Build codimension-one topology from all cells, including mixed cell types.
nCells = grid.cells.num;
halfFaces = cell(1,numel(grid.cells.vtkTypes));
halfCells = halfFaces;
halfLocal = halfFaces;

for k = 1:numel(grid.cells.vtkTypes)
  vtkId = grid.cells.vtkTypes(k);
  ids = grid.getCellsByVTKId(vtkId);
  nodes = grid.getCellNodes(ids);
  local = localFaceDefs(vtkId);
  nLocal = size(local,1);
  width = size(local,2);
  halfFaces{k} = zeros(numel(ids)*nLocal,4);
  for j = 1:nLocal
    rows = (j-1)*numel(ids) + (1:numel(ids));
    conn = nodes(:,local(j,:));
    % Canonicalize cyclic connectivity without changing the polygon perimeter.
    [~,first] = min(conn,[],2);
    offsets = mod(first-1+(0:width-1),width);
    conn = conn((1:numel(ids)).' + offsets*numel(ids));
    flip = conn(:,2) > conn(:,end);
    conn(flip,2:end) = conn(flip,end:-1:2);
    halfFaces{k}(rows,1:width) = conn;
  end
  halfCells{k} = repmat(ids(:),nLocal,1);
  halfLocal{k} = repelem((1:nLocal).',numel(ids));
end

conn = vertcat(halfFaces{:});
owners = vertcat(halfCells{:});
localIds = vertcat(halfLocal{:});
[faceNodes,~,faceId] = unique(conn,'rows');
nFaces = size(faceNodes,1);
counts = accumarray(faceId,1,[nFaces,1]);
assert(all(counts <= 2), 'GReS:NonManifoldGrid', ...
  'A face belongs to more than two cells.');
[~,order] = sortrows([faceId,(1:numel(faceId)).']);
first = cumsum([1;counts(1:end-1)]);
neighbors = [owners(order(first)),zeros(nFaces,1)];
internal = counts == 2;
neighbors(internal,2) = owners(order(first(internal)+1));

f = grid.faces;
f.num = nFaces;
f.connectivity = ArrayOfArrays(faceNodes);
f.numVerts = sum(faceNodes > 0,2);
f.neighbors = neighbors;
f.isBoundary = ~internal;
if grid.nDim == 2
  p1 = grid.coordinates(faceNodes(:,1),:);
  p2 = grid.coordinates(faceNodes(:,2),:);
  tangent = p2 - p1;
  f.area = vecnorm(tangent,2,2);
  f.center = (p1 + p2)/2;
  if grid.nCoord == 2
    tangent(:,3) = 0;
  end
  f.normal = cross(tangent,grid.cells.normal(neighbors(:,1),:),2) ./ f.area;
  f.normal = f.normal(:,1:grid.nCoord);
  f.VTKType = repmat(double(VTKType.Line),nFaces,1);
else
  [f.area,f.center,f.normal] = computePolygonGeometry( ...
    grid.coordinates(f.connectivity.getData(),:),f.numVerts);
  f.VTKType = repmat(double(VTKType.Tri),nFaces,1);
  f.VTKType(f.numVerts == 4) = double(VTKType.Quad);
end
grid.faces = f;

[~,order] = sortrows([owners,localIds]);
nLocal = accumarray(owners,1,[nCells,1]);
grid.cells.cells2faces = ArrayOfArrays(faceId(order),nLocal);
grid.cells.cells2localFaces = ArrayOfArrays(localIds(order),nLocal);

s = grid.surfaces;
boundaryIds = find(f.isBoundary);
if s.num == 0
  s.faceId = boundaryIds;
  s.connectivity = getRows(f.connectivity,boundaryIds);
  s.num = numel(boundaryIds);
  s.numVerts = f.numVerts(boundaryIds);
  s.VTKType = f.VTKType(boundaryIds);
  s.tag = ones(s.num,1);
else
  keys = zeros(s.num,4);
  for vtkId = s.vtkTypes
    ids = grid.getSurfByVTKId(vtkId);
    nodes = grid.getSurfNodes(ids);
    nodes = nodes(:,1:s.numVerts(ids(1)));
    if vtkId == VTKType.Quad9
      nodes = nodes(:,1:4);
    end
    keys(ids,1:size(nodes,2)) = sort(nodes,2);
  end
  % Map each supplied boundary row separately, including repeated tags.
  [found,face] = ismember(sort(keys,2),sort(faceNodes(boundaryIds,:),2),'rows');
  assert(all(found), 'GReS:BoundaryFaceMismatch', ...
    'Every tagged boundary entity must match a boundary face.');
  s.faceId = boundaryIds(face);
end
s.nTag = numel(unique(s.tag));
s.vtkTypes = reshape(unique(s.VTKType),1,[]);
grid.surfaces = s;
end

function faces = localFaceDefs(vtkId)
switch double(vtkId)
  case 5
    faces = [1 2; 2 3; 3 1];
  case {9,28}
    faces = [1 2; 2 3; 3 4; 4 1];
  case 10
    faces = [1 2 4; 2 3 4; 3 1 4; 1 3 2];
  case {12,29}
    faces = [1 4 3 2; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
  otherwise
    error('GReS:UnsupportedCellType','Unsupported VTK cell type %d.',vtkId);
end
end
