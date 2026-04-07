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
%
%  See also: processGeometry, ArrayOfArrays


% 1.  Extract typed half-faces for every cell type
%
%     triHF  (nTriHF  × nVert + 2) : [n1 n2 n3    cellId localFaceId]
%     quadHF (nQuadHF × nVers + 2) : [n1 n2 n3 n4 cellId localFaceId]


if nargin == 1
  % if vtk id is not provided, recursively process all vtk type available
  vtkIds = unique(grid.cells.VTKType');

  for vtkId = vtkIds
    processFaces(grid,vtkId);
  end

  return
end

coords = grid.coordinates;       
nCells = grid.cells.num;

% 1) half-face connectivity
idC = find(grid.cells.VTKType == vtkId);
topol = grid.cells.getCellNodes(idC);

f = localFaceDefs(vtkId);
nNPF = size(f,2);

% half-faces
hf = zeros(size(f,1)*nCells,size(f,2));

% repeated half-face topology
for i = 1:size(f,1)
  k = (i-1)*nCells;
  hf(k+1:k+nCells) = topol(:,f(i,:));
end

% cell index per half face
hfCellId   = repelem(idC,nNPF);

% local face id per face
hfLocalFId = repmat((1:nNPF)', nC, 1);  



% partially inspired by MRST routine

% canonicalize faces
nHF = size(hf, 1);

% -- Step 1: rotate minimum to column 1  ---------------------------------
[~, minCol] = min(hfNodes, [], 2);                             % nHF × 1, 1-based
offsets     = mod(bsxfun(@plus, (0:nNPF-1), minCol-1), nNPF); % nHF × nNPF, 0-based column offsets

% Column-major linear index: element at (row i, col offsets(i,j)+1)
% linear = offsets(i,j)*nHF + i
linIdx = bsxfun(@plus, offsets * nHF, (1:nHF)');              % nHF × nNPF
hfCan  = reshape(hfNodes(linIdx), nHF, nNPF);

% -- Step 2: choose canonical cyclic direction  --------------------------
needFlip               = hfCan(:,2) > hfCan(:,end);
hfCan(needFlip, 2:end) = hfCan(needFlip, end:-1:2);

hfSign           = ones(nHF, 1);
hfSign(needFlip) = -1;

% --  Step B: sort to bring matching half-faces to adjacent rows  --------
[hfSorted, sIdx] = sortrows(hfCan);
hfCellSorted     = hfCellId  (sIdx);
hfSignSorted     = hfSign    (sIdx);
hfLocalSorted    = hfLocalFId(sIdx);

% --  Step C: identify unique faces  -------------------------------------
%  Consecutive equal rows → same face.  First occurrence marks a new face.
isNew    = [true; any(diff(hfSorted), 2)];   % nHF × 1  logical
faceIdx  = cumsum(isNew);                     % half-face → face ID (1-based)
nF       = faceIdx(end);                     % number of faces
faceNodes = hfSorted(isNew, :);              % nF × nNPF

% --  Step D: build neighbor table  --------------------------------------
%  hfSign = +1  →  column 1 of neighbors (n1)
%  hfSign = -1  →  column 2 of neighbors (n2)
col       = ones(nHF, 1);
col(hfSignSorted == -1) = 2;
neighbors = accumarray([faceIdx, col], hfCellSorted, [nF, 2]);

% guarantee that boundary faces always have their cell in column 1
needSwap              = neighbors(:,1) == 0 & neighbors(:,2) ~= 0;
neighbors(needSwap,:) = neighbors(needSwap, [2,1]);

% --  Step E: cell-to-face map  ------------------------------------------
cellFaces = [hfCellSorted, faceIdx, hfLocalSorted];


% process face geometry 

[a,c,n] = computePolygonGeometry(poly,nNPF*ones(nF));


% update grid structure
nVerts = [3*ones(nTri,1); 4*ones(nQuad,1)];

f.num       = nF;
f.neighbors = [triNbr;  quadNbr];
f.normals   = [triNrm;  quadNrm];
f.centroids = [triCen;  quadCen];
f.areas     = [triA;    quadA  ];
f.numVerts  = nVerts;

% concatenate connectivity maps into ArrayOfArrays
if nTri > 0 && nQuad > 0
  % Mixed face types: pack into an ArrayOfArrays
  flatNodes      = [reshape(triFN',  [], 1);    % 3·nTri  entries
    reshape(quadFN', [], 1)];   % 4·nQuad entries
  f.connectivity = ArrayOfArrays(flatNodes, nVerts);
elseif nTri > 0
  f.connectivity = int32(triFN);
else
  f.connectivity = int32(quadFN);
end

grid.faces = f;

% -----------------------------------------------------------------------
% 5.  Update grid.cells cell-to-face mapping
%     allCF rows:  [cellId  globalFaceId  localFaceId]
% -----------------------------------------------------------------------
if ~isempty(quadCF)
  quadCF(:,2) = quadCF(:,2) + nTri;   % shift to global face numbering
end

allCF  = sortrows([triCF; quadCF], 1);            % sort by cell ID
nFPCv  = accumarray(allCF(:,1), 1, [nCells, 1]);  % faces per cell

grid.cells.facePos     = int32([1; cumsum(nFPCv) + 1]);
grid.cells.faces       = int32(allCF(:, 2));
grid.cells.faceLocalId = int32(allCF(:, 3));

end


% fix normals 


% map external surfaces to faces







%-------------------------------------------------------------------------

function lf = localFaceDefs(vtkId)
%LOCALFACEDEFS  Return local face node indices (1-based) for a VTK cell.
%
% Node ordering of faces is CCW when viewed from outside.

switch vtkId

  case 10
    lf   = [ 1, 2, 4 ;
      2, 3, 4 ;
      3, 1, 4 ;
      1, 3, 2 ];
  case {12, 29}
    lf   = [ 1, 4, 3, 2 ;
      5, 6, 7, 8 ;
      1, 2, 6, 5 ;
      2, 3, 7, 6 ;
      3, 4, 8, 7 ;
      4, 1, 5, 8 ];
end
end

