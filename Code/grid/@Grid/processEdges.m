function processEdges(grid,vtkId)
%PROCESSEDGES  Build edge topology and geometry for a manifold surface mesh.
%
%  PROCESSEDGES(grid) operates on all 2-D cell types found in grid.surfaces
%  and populates:
%
%  grid.edges
%    .num          – scalar, total number of unique edges
%    .connectivity – int32 matrix  (nE × 2), global node ids of each edge
%    .neighbors    – nE × 2  [s1, s2] global surface IDs; s2 = 0 on boundary
%    .length       – nE × 1
%    .center       – nE × nDim  edge midpoint
%    .isBoundary   – nE × 1  % logical true if edge is on the boundary
%
%  grid.surfaces (appended)
%    .surfaces2edges      – packed global edge IDs  (one entry per surface-edge pair)
%    .surfaces2localEdges – packed local edge index inside the surface  (1-based)
%
%  Edge convention
%    connectivity is stored in canonical order [min(nodeId) max(nodeId)]
%
%  Neighbor convention
%    boundary : neighbors(:,2) = 0
%    internal : edge shared by exactly two surfaces
%
%  Supported VTK types: triangle, quadrilateral, biquadratic quadrilateral
%  IMPORTANT: This function assumes a manifold surface mesh, so each edge
%  can belong to at most two surfaces.
%
%  See also: processFaces, ArrayOfArrays

if nargin == 1
  for vtkId = grid.surfaces.vtkTypes
    processEdges(grid, vtkId);
  end
  return
end

coords = grid.coordinates;

idS = grid.getSurfByVTKId(vtkId);
if isempty(idS)
  return
end

nS    = numel(idS);
topol = grid.getSurfNodes(idS);

le = localEdgeDefs(vtkId);
[nEPS, ~] = size(le);

% half-edges
heNodes = zeros(nEPS*nS, 2, class(topol));
for id = 1:nEPS
  k = (id-1)*nS;
  heNodes(k+1:k+nS, :) = topol(:, le(id,:));
end

heSurfId   = repmat(idS(:), nEPS, 1);
heLocalEId = repelem((1:nEPS)', nS);

% canonical representation of edges
heCan = sort(heNodes, 2);

[heSorted, sIdx] = sortrows(heCan);
heSurfSorted     = heSurfId(sIdx);
heLocalSorted    = heLocalEId(sIdx);

isNew   = [true; any(heSorted(2:end,:) ~= heSorted(1:end-1,:), 2)];
edgeIdx = cumsum(isNew);
nE      = edgeIdx(end);

% manifold check
nSPE = accumarray(edgeIdx, 1, [nE, 1]); % number of surfaces per edge
if any(nSPE > 2)
  error('processEdges:nonManifoldMesh', ...
    'Non-manifold surface mesh detected: some edges belong to more than two surfaces.');
end

edgeTopol = heSorted(isNew,:);

% get unique neighbors
col = ones(nEPS*nS, 1);
col(~isNew) = 2;
neighbors = accumarray([edgeIdx, col], heSurfSorted, [nE, 2]);

% make sure that boundary edges have 0 surface as second position
needSwap = neighbors(:,1) == 0 & neighbors(:,2) ~= 0;
neighbors(needSwap,:) = neighbors(needSwap, [2,1]);

isBoundary = neighbors(:,2) == 0;

% compute edge geometrical informations
p1 = coords(edgeTopol(:,1),:);
p2 = coords(edgeTopol(:,2),:);

center = 0.5*(p1 + p2);
len = sqrt(sum((p2 - p1).^2, 2));

% append edges
e = grid.edges;
s = grid.surfaces;

nEold = e.num;
e.num          = nEold + nE;
e.connectivity = [e.connectivity; edgeTopol];
e.neighbors    = [e.neighbors; neighbors];
e.length       = [e.length; len];
e.center       = [e.center; center];
e.isBoundary   = [e.isBoundary; isBoundary];

% surface to edges mapping
[~,id] = sort(heSurfSorted);
nEPSs = nEPS*ones(nS,1);

s2e = ArrayOfArrays(nEold + edgeIdx(id), nEPSs);
%s2loce = ArrayOfArrays(heLocalSorted(id), nEPSs);

s.surfaces2edges = [s.surfaces2edges; s2e];
%s.surfaces2localEdges = [s.surfaces2localEdges; s2loce];

% finally update grid
grid.edges = e;
grid.surfaces = s;

end


function le = localEdgeDefs(vtkId)

switch double(vtkId)
  case 5
    le = [1 2;
          2 3;
          3 1];

  case {9, 28}
    le = [1 2;
          2 3;
          3 4;
          4 1];

  otherwise
    error('processEdges:unsupportedVTKType', ...
      'Unsupported VTK type %d.', vtkId);
end

end