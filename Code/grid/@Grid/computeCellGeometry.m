function [volumes, centers] = computeCellGeometry(grid)
%COMPUTECELLGEOMETRY Compute cell measures and centroids in 2D or 3D
%
%   [volumes, centers] = computeCellGeometry(grid)
%
% Geometry is computed with the divergence theorem.
% volumes contains areas for grid.nDim = 2 and volumes for grid.nDim = 3.
% Coordinates and returned centers may have two or three components in 2D.
%
% Assumptions:
%   - grid.nDim                     : intrinsic dimension, 2 or 3
%   - grid.coordinates              : nN x nCoord
%   - grid.cells.num                : number of cells
%   - grid.cells.connectivity       : matrix or ArrayOfArrays-compatible
%   - grid.cells.cells2faces        : matrix or ArrayOfArrays-compatible
%   - grid.faces.connectivity      : edge endpoint nodes in 2D
%   - grid.faces.center             : face centroids or edge midpoints
%   - grid.faces.normal             : unit face normals, used in 3D only
%   - grid.faces.area               : face areas in 3D, edge lengths in 2D
%
% Notes:
%   - In 2D, each cell must be planar and its edges straight.
%   - The rough centroid must lie on the interior side of every supporting
%     face/edge. Convex cells satisfy this condition.
%   - In 2D, distances are computed separately for each cell-edge pair,
%     allowing neighboring cells to lie in different planes.
%   - The returned centers are the final geometric centroids.

coords = grid.coordinates;
nC = grid.cells.num;
nDim = grid.nDim;
nCoord = size(coords,2);

cellConn = ArrayOfArrays(grid.cells.connectivity);
[cellNodes,ptr] = getData(cellConn);
nNPC = diff(ptr);

xyz = coords(cellNodes,:);
rowId = repelem((1:nC).', nNPC, 1);

roughCenters = zeros(nC, nCoord);
for k = 1:nCoord
  roughCenters(:,k) = accumarray(rowId, xyz(:,k), [nC,1]) ./ nNPC;
end


cell2faces = ArrayOfArrays(grid.cells.cells2faces);
[c2f, ptr] = getData(cell2faces);
nFPC = diff(ptr);

% expand cells for each face
rowId = repelem((1:nC).', nFPC, 1);


cf = grid.faces.center(c2f,:);
cc = roughCenters(rowId,:);
relcf = cf - cc;

if nDim == 2
  % Distance to each supporting edge line within the cell plane.
  % This avoids using one shared co-normal for noncoplanar cells.
  faceConn = ArrayOfArrays(grid.faces.connectivity);
  [faceNodes,facePtr] = getData(faceConn);

  p1 = coords(faceNodes(facePtr(c2f)),:);
  p2 = coords(faceNodes(facePtr(c2f) + 1),:);
  tangent = p2 - p1;
  tangent = tangent ./ sqrt(sum(tangent.^2,2));

  normalPart = relcf - sum(relcf .* tangent,2) .* tangent;
  contribVol = sqrt(sum(normalPart.^2,2)) .* grid.faces.area(c2f);
else
  % area-weighted normals
  Nf = grid.faces.normal(c2f,:) .* grid.faces.area(c2f);

  % local signed contribution relative to rough interior point
  projLocal = sum(relcf .* Nf, 2);

  % use abs assuming cc is on the interior side of each face
  contribVol = abs(projLocal);
end

volumes = accumarray(rowId, contribVol, [nC,1]) / nDim;

% For consistency with the same local origin cc:
%   V = 1/nDim sum_f contribVol_f
%   c = cc + 1/((nDim+1)*V) sum_f contribVol_f (cf-cc)

% Compute centroid location relative to rough centroids
mom = relcf .* contribVol;

centers = zeros(nC, nCoord);
for k = 1:nCoord
  centers(:,k) = accumarray(rowId, mom(:,k), [nC,1]);
end

centers = roughCenters + centers ./ ((nDim + 1) * volumes);

end
