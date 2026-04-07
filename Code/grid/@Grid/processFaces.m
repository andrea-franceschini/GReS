function processFaces(grid,vtkId)
%PROCESSFACES  Build face topology and geometry for an unstructured grid.
%
%  PROCESSFACES(grid) operates on all 3‑D cell types found in grid and
%  populates the following fields in-place:
%
%  grid.faces
%    .num          – scalar, total number of unique faces
%    .connectivity – int32 matrix  (nF × nNPF)  for uniform meshes, or
%                    ArrayOfArrays              for mixed tet+hex meshes
%    .numVerts     – nF × 1  nodes per face  (3 = triangle, 4 = quad)
%    .neighbors    – nF × 2  [n1, n2] global cell IDs;  n2 = 0 on boundary
%    .normals      – nF × 3  area‑weighted face normals  (‖n‖₂ = area)
%    .centroids    – nF × 3  face centroid coordinates
%    .areas        – nF × 1
%    .isBoundary   – nF × 1  % logical true if face is at the boundary
%
%  grid.cells (appended)
%    .facePos      – (nCells+1) × 1  offsets into cells.faces (1‑based)
%    .faces        – packed global face IDs  (one entry per cell–face pair)
%    .faceLocalId  – packed local face index inside the cell  (1‑based)
% 
%  grid.surfaces (appended)
%    .faceId      – map tagged surface to global id of corresponding face

%  Normal convention (normal is unit)
%    boundary : normal points outward  from neighbors(:,1)
%    internal : normal points          from neighbors(:,1) to neighbors(:,2)

% Supported VTK types: tetrahedron and hexahedron (hexa27 are treated as hexa8)
% IMPORTANT: The code assumes that the processed vtkId has constant number of face (not valid for VTK type POLYHEDRON)
%  See also: processGeometry, ArrayOfArrays


% 1.  Extract typed half-faces for every cell type
%
%     triHF  (nTriHF  × nVert + 2) : [n1 n2 n3    cellId localFaceId]
%     quadHF (nQuadHF × nVers + 2) : [n1 n2 n3 n4 cellId localFaceId]


if nargin == 1
  vtkIds = unique(grid.cells.VTKType);
  for k = 1:numel(vtkIds)
    processFaces(grid, vtkIds(k));
  end
  return
end

coords = grid.coordinates; 

idC = find(grid.cells.VTKType == vtkId);
if isempty(idC)
  return
end

nC    = numel(idC);
topol = grid.getCellNodes(idC);

lf = localFaceDefs(vtkId);
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

[~, minCol] = min(hfNodes, [], 2);
offsets     = mod(bsxfun(@plus, 0:nNPF-1, minCol-1), nNPF);
linIdx      = bsxfun(@plus, offsets*nHF, (1:nHF)');
hfCan       = reshape(hfNodes(linIdx), nHF, nNPF);

needFlip               = hfCan(:,2) > hfCan(:,end);
hfCan(needFlip, 2:end) = hfCan(needFlip, end:-1:2);

% winding sign
hfSign           = ones(nHF, 1);
hfSign(needFlip) = -1;

[hfSorted, sIdx] = sortrows(hfCan);
hfCellSorted     = hfCellId(sIdx);
hfSignSorted     = hfSign(sIdx);
hfLocalSorted    = hfLocalFId(sIdx);

isNew   = [true; any(hfSorted(2:end,:) ~= hfSorted(1:end-1,:), 2)];
faceIdx = cumsum(isNew);
nF      = faceIdx(end);
nPlist = nNPF*ones(nF,1);

faceTopol = int32(hfSorted(isNew, :));

% get unique neighbors
col = ones(nHF, 1);
col(hfSignSorted == -1) = 2;
neighbors = accumarray([faceIdx, col], hfCellSorted, [nF, 2]);

% make sure that boundary faces have 0 cell as secondo position
needSwap = neighbors(:,1) == 0 & neighbors(:,2) ~= 0;
neighbors(needSwap,:) = neighbors(needSwap, [2,1]);

isBoundary = neighbors(:,2) == 0;

% compute faces geometrical informations
nList = faceTopol';
poly = coords(nList(:),:); 
[area, cent, normal] = computePolygonGeometry(poly, nPlist);


% append faces
f = grid.faces;
c = grid.cells;
s = grid.surfaces;

nFold = f.num;
f.num        = nFold + nF;
f.neighbors  = [f.neighbors; int32(neighbors)];
f.normal    = [f.normal; normal];
f.center  = [f.center; cent];
f.area      = [f.area; area];
f.numVerts   = [f.numVerts; nPlist];
f.isBoundary = [f.isBoundary; isBoundary];
f.connectivity = [f.connectivity; ArrayOfArrays(faceTopol)];


% cell to faces mapping
[~,id] = sort(hfCellSorted);
nFPCs = nFPC*ones(nC,1);
c2f = ArrayOfArrays(faceIdx(id),nFPCs);
c2locf = ArrayOfArrays(hfLocalSorted(id),nFPCs);
c.cells2faces = [c.cells2faces; c2f];
c.cells2localFaces = [c.cells2localFaces; c2locf];

% surfaces (map surfaces to boundaryfaces) 
boundFaces = faceTopol(isBoundary,:);

% connectivity of the surfaces corresponding to the 2D VTK
vtk2d = grid.vtkType(grid.vtkType(:,2) == vtkId,1);
idS = find(grid.surfaces.VTKType == vtk2d);
topolSurf = grid.getSurfNodes(idS);
[sId,fId] = ismember(sort(topolSurf,2),sort(boundFaces,2),'rows');
boundId = find(isBoundary);
s.faceId(sId) = nFold + boundId(fId); 


% finally update grid
grid.faces = f;
grid.cells = c;
grid.surfaces = s;


end


function lf = localFaceDefs(vtkId)

switch vtkId
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

