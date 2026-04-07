function [volumes, centers] = computeCellGeometry(grid)
%COMPUTECELLGEOMETRY Compute cell volumes and centroids for generic polyhedra
%
%   [volumes, centers] = computeCellGeometry(grid)
%
% Geometry is computed with the divergence theorem. 
% It is assumed that the average of the vertex coordinates (rough centroid)
% lies inside the cell
%
% Assumptions:
%   - grid.coordinates              : nN x 3
%   - grid.cells.num                : number of cells
%   - grid.cells.connectivity       : matrix or ArrayOfArrays-compatible
%   - grid.cells.cell2faces         : matrix or ArrayOfArrays-compatible
%   - grid.faces.center             : nF x 3
%   - grid.faces.normal             : nF x 3  (unit normals)
%   - grid.faces.area               : nF x 1
%
% Notes:
%   - This computes geometric polyhedral volume/centroid.
%   - The rough centroid is only used to define local outwardness through abs().
%   - The returned centers are the final geometric centroids, not the rough ones.

coords = grid.coordinates;
nC = grid.cells.num;

cellConn = ArrayOfArrays(grid.cells.connectivity);
[cellNodes,ptr] = getData(cellConn);   % flat node ids, nodes-per-cell
nNPC = diff(ptr);

xyz = coords(cellNodes,:);
rowId = repelem((1:nC).', nNPC);

roughCenters = zeros(nC, 3);
roughCenters(:,1) = accumarray(rowId, xyz(:,1), [nC,1]) ./ nNPC;
roughCenters(:,2) = accumarray(rowId, xyz(:,2), [nC,1]) ./ nNPC;
roughCenters(:,3) = accumarray(rowId, xyz(:,3), [nC,1]) ./ nNPC;


cell2faces = ArrayOfArrays(grid.cells.cells2faces);
[c2f, ptr] = getData(cell2faces);    
nFPC = diff(ptr);

% expand cells for each face
rowId = repelem((1:nC).', nFPC);


cf = grid.faces.center(c2f,:);           % face centroids

 % area-weighted normals
Nf = grid.faces.normal(c2f,:) .* grid.faces.area(c2f);  

cc = roughCenters(rowId,:);


% local signed contribution relative to rough interior point
projLocal = sum((cf - cc) .* Nf, 2);

% use abs assuming cc lies within the polyhedron
contribVol = abs(projLocal);

volumes = accumarray(rowId, contribVol, [nC,1]) / 3;

% For consistency with the same local origin cc:
%   V  = 1/3 sum_f |(cf-cc)·Nf|
%   c  = cc + (1/(4V)) sum_f |(cf-cc)·Nf| (cf-cc)
%

% Comptute centroid location relative to rough centroids

relcf = cf - cc;
mom = relcf .* contribVol;   % each row scaled by its scalar contribution

Mx = accumarray(rowId, mom(:,1), [nC,1]);
My = accumarray(rowId, mom(:,2), [nC,1]);
Mz = accumarray(rowId, mom(:,3), [nC,1]);

centers = roughCenters + [Mx, My, Mz] ./ (4 * volumes);

end