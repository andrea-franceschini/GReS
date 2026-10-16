function processFaces(grid,vtkId,varargin)
%PROCESSFACES Build face or edge topology and geometry.
%
% processFaces(grid,vtkId)
%   nDim = 3: process grid.cells into grid.faces (polygonal faces).
%   nDim = 2: process grid.cells into grid.faces (edges).
%
% processFaces(grid,vtkId,"edges")
%   Process grid.surfaces into grid.edges. Requires grid.nDim = 3.
%
% If vtkId is omitted or empty, loop over the source VTK types.
% Each VTK type is processed separately
% and its topology and geometry are appended to the existing fields.
%
% In 2D face mode, face.area is edge length and face.normal is a unit
% co-normal computed in the plane of neighbors(:,1). Cells must be planar.
% Higher-order elements are treated using their corner nodes.

nDim = grid.nDim;
isEdges = ~isempty(varargin) && strcmpi(varargin{1},'edges');

if isEdges && nDim ~= 3
  error('processFaces:invalidDimension', ...
    'The "edges" option requires grid.nDim = 3.');
end

if isEdges
  source = grid.surfaces;
else
  source = grid.cells;
end

if nargin == 1 || isempty(vtkId)
  vtkIds = unique(source.VTKType);
  for k = 1:numel(vtkIds)
    if isEdges
      processFaces(grid, vtkIds(k), "edges");
    else
      processFaces(grid, vtkIds(k));
    end
  end
  return
end

coords = grid.coordinates;
isLine = isEdges || nDim == 2;

if isEdges
  idC = grid.getSurfByVTKId(vtkId);
else
  idC = grid.getCellsByVTKId(vtkId);
end

if isempty(idC)
  return
end

nC = numel(idC);
if isEdges
  topol = grid.getSurfNodes(idC);
else
  topol = grid.getCellNodes(idC);
end

lf = localFaceDefs(vtkId,isLine);
[nFPC, nNPF] = size(lf);

% half-faces
hfNodes = zeros(nFPC*nC, nNPF, class(topol));
for id = 1:nFPC
  k = (id-1)*nC;
  hfNodes(k+1:k+nC, :) = topol(:, lf(id,:));
end

hfCellId   = repmat(idC(:), nFPC, 1);
hfLocalFId = repelem((1:nFPC)', nC);

% canonical representation of faces
nHF = size(hfNodes, 1);

if isLine
  hfCan = sort(hfNodes, 2);
else
  [~, minCol] = min(hfNodes, [], 2);
  offsets = mod(bsxfun(@plus, 0:nNPF-1, minCol-1), nNPF);
  linIdx = bsxfun(@plus, offsets*nHF, (1:nHF)');
  hfCan = reshape(hfNodes(linIdx), nHF, nNPF);

  needFlip = hfCan(:,2) > hfCan(:,end);
  hfCan(needFlip, 2:end) = hfCan(needFlip, end:-1:2);
end

[hfSorted, sIdx] = sortrows(hfCan);
hfCellSorted = hfCellId(sIdx);
hfLocalSorted = hfLocalFId(sIdx);

isNew = [true; any(hfSorted(2:end,:) ~= hfSorted(1:end-1,:), 2)];
faceIdx = cumsum(isNew);
nF = faceIdx(end);
nPlist = nNPF*ones(nF,1);

faceTopol = hfSorted(isNew, :);

% get unique neighbors
nCPF = accumarray(faceIdx, 1, [nF, 1]);
if any(nCPF > 2)
  error('processFaces:nonManifoldMesh', ...
    'An entity belongs to more than two source elements.');
end

col = ones(nHF, 1);
col(~isNew) = 2;
neighbors = accumarray([faceIdx, col], hfCellSorted, [nF, 2]);

% make sure that boundary faces have 0 cell as second position
needSwap = neighbors(:,1) == 0 & neighbors(:,2) ~= 0;
neighbors(needSwap,:) = neighbors(needSwap, [2,1]);

isBoundary = neighbors(:,2) == 0;

% compute faces geometrical informations
if isLine
  p1 = coords(faceTopol(:,1),:);
  p2 = coords(faceTopol(:,2),:);

  cent = 0.5*(p1 + p2);
  area = sqrt(sum((p2 - p1).^2, 2));

  if ~isEdges
    % Unit co-normal in the plane of the first neighboring 2D cell.
    % Surface corner nodes are in cyclic order, including for Quad9.
    nCorners = nFPC;
    cellTopol = topol(:,1:nCorners);
    nList = cellTopol';
    poly = coords(nList(:),:);
    if size(poly,2) == 2
      poly(:,3) = 0;
    end
    [~, cellCenter, cellNormal] = computePolygonGeometry(poly, nCorners*ones(nC,1));

    [~, id] = ismember(neighbors(:,1),idC);
    tangent = (p2 - p1) ./ area;
    if size(tangent,2) == 2
      tangent(:,3) = 0;
    end
    normal = cross(tangent,cellNormal(id,:),2);

    direction = zeros(nF,3);
    direction(:,1:size(coords,2)) = cent;
    direction = direction - cellCenter(id,:);
    needFlip = sum(normal .* direction,2) < 0;
    normal(needFlip,:) = -normal(needFlip,:);
    normal = normal(:,1:size(coords,2));
  end
else
  nList = faceTopol';
  poly = coords(nList(:),:);
  [area, cent, normal] = computePolygonGeometry(poly, nPlist);
end

% append faces or edges
if isEdges
  f = grid.edges;
else
  f = grid.faces;
end

nFold = f.num;
f.num = nFold + nF;
f.neighbors = [f.neighbors; neighbors];
f.center = [f.center; cent];
f.isBoundary = [f.isBoundary; isBoundary];

if isEdges
  f.connectivity = [f.connectivity; faceTopol];
  f.length = [f.length; area];
else
  f.normal = [f.normal; normal];
  f.area = [f.area; area];
  f.numVerts = [f.numVerts; nPlist];
  f.connectivity = [f.connectivity; ArrayOfArrays(faceTopol)];
  if nDim == 3
    f.VTKType = double(VTKType.to2D(source.VTKType(f.neighbors(:,1))));
  else
    f.VTKType = 3*ones(f.num,1);
  end
end

% source element to faces mapping
[~,id] = sort(hfCellSorted);
nFPCs = nFPC*ones(nC,1);
c2f = ArrayOfArrays(nFold + faceIdx(id),nFPCs);
c2locf = ArrayOfArrays(hfLocalSorted(id),nFPCs);

if isEdges
  source.surfaces2edges = [source.surfaces2edges; c2f];
  if isfield(source,'surfaces2localEdges')
    source.surfaces2localEdges = [source.surfaces2localEdges; c2locf];
  else
    source.surfaces2localEdges = c2locf;
  end
else
  source.cells2faces = [source.cells2faces; c2f];
  source.cells2localFaces = [source.cells2localFaces; c2locf];
end

% surfaces (map surfaces to boundary faces)
if ~isEdges
  s = grid.surfaces;
  boundFaces = faceTopol(isBoundary,:);

  if s.num > 0
    if nDim == 3
      surfVtkId = VTKType.to2D(vtkId);
    else
      surfVtkId = 3;   % VTK_LINE
    end
    idS = grid.getSurfByVTKId(surfVtkId);
    if ~isempty(idS)
      topolSurf = grid.getSurfNodes(idS);
      if nDim == 3 && surfVtkId == VTKType.Quad9
        topolSurf = topolSurf(:,1:4);
      end
      [sId,fId] = ismember(sort(topolSurf,2),sort(boundFaces,2),'rows');
      boundId = find(isBoundary);
      s.faceId(idS(sId)) = nFold + boundId(fId(sId));
    end
  end

  grid.surfaces = s;
end


% finally update grid
if isEdges
  grid.edges = f;
  grid.surfaces = source;
else
  grid.faces = f;
  grid.cells = source;
end

end


function lf = localFaceDefs(vtkId,isLine)

if isLine
  switch double(vtkId)
    case 5
      lf = [1 2;
        2 3;
        3 1];

    case {9, 28}
      lf = [1 2;
        2 3;
        3 4;
        4 1];

    otherwise
      error('processFaces:unsupportedVTKType', ...
        'Unsupported VTK type %d.', vtkId);
  end
else
  switch double(vtkId)
    case 10
      lf = [1 2 4;
        2 3 4;
        3 1 4;
        1 3 2];

    case {12, 29}
      lf = [1 4 3 2;
        5 6 7 8;
        1 2 6 5;
        2 3 7 6;
        3 4 8 7;
        4 1 5 8];

    otherwise
      error('processFaces:unsupportedVTKType', ...
        'Unsupported VTK type %d.', vtkId);
  end
end

end
