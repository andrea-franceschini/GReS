function processFaces(grid)
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
%
%  grid.cells (appended)
%    .facePos      – (nCells+1) × 1  offsets into cells.faces (1‑based)
%    .faces        – packed global face IDs  (one entry per cell–face pair)
%    .faceLocalId  – packed local face index inside the cell  (1‑based)
%
%  Normal convention (normal is unit)
%    boundary : normal points outward  from neighbors(:,1)
%    internal : normal points          from neighbors(:,1) to neighbors(:,2)

%
%  Supported VTK 3‑D cell types
%    10   VTK_TETRA                     4 nodes  → 4 triangular faces
%    12   VTK_HEXAHEDRON                8 nodes  → 6 quad faces
%    29   VTK_TRIQUADRATIC_HEXAHEDRON  27 nodes  → 6 quad faces (corner nodes)
%
%  See also: processGeometry, ArrayOfArrays

  coords = double(grid.coordinates);         % nNodes × 3
  nCells = numel(grid.cells.VTKType);

  % -----------------------------------------------------------------------
  % 1.  Extract typed half-faces for every cell type
  %
  %     triHF  (nTriHF  × 5) : [n1 n2 n3    cellId localFaceId]
  %     quadHF (nQuadHF × 6) : [n1 n2 n3 n4 cellId localFaceId]
  %
  %     A "half-face" is a face as seen from one cell.  Two half-faces that
  %     share the same node set (in opposite winding) form one unique face.
  % -----------------------------------------------------------------------
  triHF  = zeros(0, 5);
  quadHF = zeros(0, 6);

  for vtkId = reshape(unique(double(grid.cells.VTKType)), 1, [])

    mask    = double(grid.cells.VTKType) == vtkId;
    cellIds = find(mask);                            % global cell IDs
    conn    = double(grid.getCellNodes(cellIds));    % nC × nNodesPerCell
    nC      = numel(cellIds);

    [lf, nNPF] = localFaceDefs(vtkId);              % lf : nFPC × nNPF
    nFPC       = size(lf, 1);                        % faces per cell

    % --  Vectorised half-face extraction  --
    %
    % We want:
    %   hfNodes((c-1)*nFPC + f, k)  =  conn(c, lf(f, k))
    %
    % For each node position k:
    %   tmp = conn(:, lf(:,k))  is  nC × nFPC,  tmp(c,f) = conn(c, lf(f,k))
    %   reshape(tmp', [], 1)  reads tmp' column-major  ≡  tmp row-major:
    %   result position (c-1)*nFPC + f  =  conn(c, lf(f,k))   ✓
    hfNodes = zeros(nC * nFPC, nNPF);
    for k = 1 : nNPF
      tmp          = conn(:, lf(:, k));          % nC × nFPC
      hfNodes(:,k) = reshape(tmp', [], 1);
    end

    hfCellId   = repelem(cellIds(:),  nFPC);     % (nC·nFPC) × 1
    hfLocalFId = repmat((1:nFPC)', nC, 1);       % local face index 1…nFPC, repeated per cell

    if nNPF == 3
      triHF  = [triHF;  hfNodes, hfCellId, hfLocalFId]; %#ok<AGROW>
    else
      quadHF = [quadHF; hfNodes, hfCellId, hfLocalFId]; %#ok<AGROW>
    end
  end

  % -----------------------------------------------------------------------
  % 2.  Build unique faces, neighbors and geometry — one pass per shape
  % -----------------------------------------------------------------------
  [triFN,  triNbr,  triNrm,  triCen,  triA,  triCF]  = ...
      buildGroup(triHF,  3, coords, nCells);
  nTri = size(triFN, 1);

  [quadFN, quadNbr, quadNrm, quadCen, quadA, quadCF] = ...
      buildGroup(quadHF, 4, coords, nCells);
  nQuad = size(quadFN, 1);
  nF    = nTri + nQuad;

  % -----------------------------------------------------------------------
  % 3.  Orient normals: outward for boundary / n1→n2 for internal
  %     Uses rough cell centroids (mean of corner nodes) for orientation.
  % -----------------------------------------------------------------------
  cellCC = roughCellCentroids(grid, coords);

  [triNrm,  triNbr]  = orientNormals(triNrm,  triCen,  triNbr,  cellCC);
  [quadNrm, quadNbr] = orientNormals(quadNrm, quadCen, quadNbr, cellCC);

  % -----------------------------------------------------------------------
  % 4.  Assemble grid.faces
  % -----------------------------------------------------------------------
  nVerts = [3*ones(nTri,1); 4*ones(nQuad,1)];

  f.num       = nF;
  f.neighbors = [triNbr;  quadNbr];
  f.normals   = [triNrm;  quadNrm];
  f.centroids = [triCen;  quadCen];
  f.areas     = [triA;    quadA  ];
  f.numVerts  = nVerts;

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


% =========================================================================
%  buildGroup
%  Core routine: deduplicates typed half-faces, infers neighbors,
%  and computes geometry for one face shape (tri or quad).
% =========================================================================
function [faceNodes, neighbors, normals, centroids, areas, cellFaces] = ...
    buildGroup(hf, nNPF, coords, nCells)
%BUILDGROUP  Deduplicate half-faces, build neighbor table, compute geometry.
%
%  hf       : nHF × (nNPF+2)  — columns: [node1…nodeN | cellId | localFaceId]
%  nNPF     : nodes per face   (3 = triangle, 4 = quad)
%  coords   : nNodes × 3
%  nCells   : total cell count (upper bound for accumarray)
%
%  faceNodes : nF × nNPF   canonical face node connectivity
%  neighbors : nF × 2      [n1, n2]  global cell IDs; 0 on boundary
%  cellFaces : nHF × 3     [cellId  globalFaceId  localFaceId]

  if isempty(hf)
    faceNodes = zeros(0, nNPF);
    neighbors = zeros(0, 2);
    normals   = zeros(0, 3);
    centroids = zeros(0, 3);
    areas     = zeros(0, 1);
    cellFaces = zeros(0, 3);
    return
  end

  hfNodes    = double(hf(:, 1:nNPF));
  hfCellId   = double(hf(:, nNPF+1));
  hfLocalFId = double(hf(:, nNPF+2));
  nHF        = size(hfNodes, 1);

  % --  Step A: canonicalize half-face node ordering  ----------------------
  %
  % Two half-faces that share a geometric face will have the same canonical
  % row but opposite winding, which we track with hfSign:
  %   +1  canonical ordering preserves the original winding →
  %       the cross-product normal of the canonical row points OUTWARD
  %       from the cell that contributed this half-face  → cell is n1
  %   -1  canonical ordering reverses  the original winding →
  %       normal points INWARD  to that cell              → cell is n2
  [hfCan, hfSign] = canonicalize(hfNodes, nNPF);

  % --  Step B: sort to bring matching half-faces to adjacent rows  --------
  [hfSorted, sIdx] = sortrows(hfCan);
  hfCellSorted     = hfCellId  (sIdx);
  hfSignSorted     = hfSign    (sIdx);
  hfLocalSorted    = hfLocalFId(sIdx);

  % --  Step C: identify unique faces  -------------------------------------
  %  Consecutive equal rows → same face.  First occurrence marks a new face.
  isNew    = [true; any(diff(hfSorted), 2)];   % nHF × 1  logical
  faceIdx  = cumsum(isNew);                     % half-face → face ID (1-based)
  nF       = faceIdx(end);
  faceNodes = hfSorted(isNew, :);              % nF × nNPF

  % --  Step D: build neighbor table  --------------------------------------
  %  hfSign = +1  →  column 1 of neighbors (n1)
  %  hfSign = -1  →  column 2 of neighbors (n2)
  col       = ones(nHF, 1);
  col(hfSignSorted == -1) = 2;
  neighbors = accumarray([faceIdx, col], hfCellSorted, [nF, 2]);

  % Guarantee that boundary faces always have their cell in column 1
  % (can occur when the sole half-face had hfSign = -1)
  needSwap              = neighbors(:,1) == 0 & neighbors(:,2) ~= 0;
  neighbors(needSwap,:) = neighbors(needSwap, [2,1]);

  % --  Step E: cell-to-face map  ------------------------------------------
  cellFaces = [hfCellSorted, faceIdx, hfLocalSorted];

  % --  Step F: geometry  --------------------------------------------------
  if nNPF == 3
    [normals, centroids, areas] = triGeom(faceNodes, coords);
  else
    [normals, centroids, areas] = quadGeom(faceNodes, coords);
  end
end


% =========================================================================
%  canonicalize
% =========================================================================
function [hfCan, hfSign] = canonicalize(hfNodes, nNPF)
%CANONICALIZE  Compute a unique representative row for each face node set.
%
%  Step 1  Cyclically rotate each row so its minimum node ID is in col 1.
%  Step 2  If col 2 > col end, reverse cols 2:end   (flip cyclic direction).
%
%  The resulting row is the same for two half-faces that share a face,
%  regardless of which cell or which starting node was used.
%
%  hfSign(i) = +1   canonical = same  winding as original
%            = -1   canonical = flipped winding  (see buildGroup for meaning)

  nHF = size(hfNodes, 1);

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
end


% =========================================================================
%  triGeom  —  vectorised geometry for triangular faces
% =========================================================================
function [normals, centroids, areas] = triGeom(faceNodes, coords)
%TRIGEOM  Area-weighted normals, centroids and areas for triangular faces.
%
%  normal(f,:) = 0.5 * cross(p2-p1, p3-p1)   →  ‖normal‖ = area
%  centroid    = (p1+p2+p3) / 3

  p1 = coords(faceNodes(:,1), :);   % nF × 3
  p2 = coords(faceNodes(:,2), :);
  p3 = coords(faceNodes(:,3), :);

  normals   = cross(p2 - p1, p3 - p1, 2) * 0.5;
  areas     = sqrt(sum(normals.^2, 2));
  centroids = (p1 + p2 + p3) / 3;

  normals = normals./areas;
end


% =========================================================================
%  quadGeom  —  vectorised geometry for (possibly non-planar) quad faces
% =========================================================================
function [normals, centroids, areas] = quadGeom(faceNodes, coords)
%QUADGEOM  Area-weighted normals, centroids and areas for quadrilateral faces.
%
%  Efficent method to compute quadrilateral geometry (assumed planar)
%  Use subtriangle subdivision
%
%  Face normal   = sum of sub-triangle normals           (area-weighted)
%  Face centroid = area-weighted mean of sub-centroids
%  Face area     = sum of sub-triangle areas

  nF   = size(faceNodes, 1);
  p    = cell(1, 4);
  for k = 1 : 4
    p{k} = coords(faceNodes(:,k), :);   % nF × 3
  end
  fCen = (p{1} + p{2} + p{3} + p{4}) * 0.25;  % coarse centroid

  next    = [2, 3, 4, 1];
  normals = zeros(nF, 3);
  areas   = zeros(nF, 1);
  wCen    = zeros(nF, 3);

  for k = 1 : 4
    sn    = cross(p{next(k)} - p{k},  fCen - p{k},  2) * 0.5;   % sub-normal
    sa    = sqrt(sum(sn.^2, 2));                                % sub-area
    sc    = (p{k} + p{next(k)} + fCen) / 3;                     % sub-centroid

    normals = normals + sn;
    areas   = areas   + sa;
    wCen    = wCen    + sa .* sc;
  end

  % guard against zero-area degenerate faces
  safeA     = max(areas, eps);   
  centroids = wCen ./ safeA;
end


% =========================================================================
%  orientNormals
% =========================================================================
function [normals, neighbors] = orientNormals(normals, centroids, neighbors, cellCC)
%ORIENTNORMALS  Enforce the outward/n1→n2 normal convention.
%
%  For internal faces the dot product of the raw normal with the vector
%  (centroid_n2 - centroid_n1) must be positive; if not, flip the normal
%  and swap the neighbor pair.
%
%  For boundary faces the dot product of the raw normal with the vector
%  (face_centroid - cell_centroid) must be positive; if not, flip.

  if isempty(normals), return; end

  isInt = neighbors(:,2) ~= 0;

  % -- Internal faces  -----------------------------------------------------
  iIdx = find(isInt);
  if ~isempty(iIdx)
    cc1  = cellCC(neighbors(iIdx,1), :);
    cc2  = cellCC(neighbors(iIdx,2), :);
    flip = sum((cc2 - cc1) .* normals(iIdx,:), 2) < 0;
    if any(flip)
      fIdx              = iIdx(flip);
      normals(fIdx,:)   = -normals(fIdx,:);
      neighbors(fIdx,:) =  neighbors(fIdx, [2,1]);
    end
  end

  % -- Boundary faces  -----------------------------------------------------
  bIdx = find(~isInt);
  if ~isempty(bIdx)
    cc   = cellCC(neighbors(bIdx,1), :);
    flip = sum((centroids(bIdx,:) - cc) .* normals(bIdx,:), 2) < 0;
    if any(flip)
      normals(bIdx(flip),:) = -normals(bIdx(flip),:);
    end
  end
end


% =========================================================================
%  roughCellCentroids
% =========================================================================
function cc = roughCellCentroids(grid, coords)
%ROUGHCELLCENTROIDS  Mean of corner-node coordinates per cell.
%
%  Higher-order elements (VTK 29) use only their 8 corner nodes so that
%  mid-edge / mid-face / body nodes do not skew the centroid estimate.
%  The centroid is only used for normal orientation, so precision is not
%  critical.

  nCells = numel(grid.cells.VTKType);
  cc     = zeros(nCells, 3);

  for vtkId = reshape(unique(double(grid.cells.VTKType)), 1, [])
    mask    = double(grid.cells.VTKType) == vtkId;
    cellIds = find(mask);
    nC      = numel(cellIds);

    conn  = double(grid.getCellNodes(cellIds));   % nC × nNodesPerCell
    nCorn = numCorners(vtkId);
    conn  = conn(:, 1:nCorn);                     % keep corner nodes only

    % Vectorised mean:
    %   reshape(conn', [], 1)  →  row-major flattening of conn
    %       position (c-1)*nCorn + k  =  conn(c, k)
    %   coords( ... )' reshaped to (3 × nCorn × nC)
    %       entry (d, k, c)  =  coords(conn(c,k), d)
    %   mean along dim 2  →  (3 × 1 × nC)  →  reshape to (3 × nC)'
    flatIdx         = reshape(conn', [], 1);            % nC·nCorn × 1
    coo             = reshape(coords(flatIdx, :)', 3, nCorn, nC);
    cc(cellIds, :)  = reshape(mean(coo, 2), 3, nC)';   % nC × 3
  end
end


% =========================================================================
%  localFaceDefs
% =========================================================================
function [lf, nNPF] = localFaceDefs(vtkId)
%LOCALFACEDEFS  Return local face node indices (1-based) for a VTK cell.
%
%  Each row of lf defines one face.  Node ordering is chosen so that the
%  cross product  (p2-p1)×(p3-p1)  points generally outward.
%  The orientNormals step will correct any remaining sign errors.
%
%  Tetrahedron  (VTK 10)  —  4 triangular faces
%  ─────────────────────
%         4
%        /|\
%       / | \
%      1--+--3
%       \ | /
%        \|/
%         2
%
%  Hexahedron  (VTK 12/29)  —  6 quad faces, corner nodes 1-8
%  ────────────────────────
%     8───────7
%    /|      /|
%   5───────6 |
%   | 4─────|─3
%   |/      |/
%   1───────2

  switch vtkId

    case 10          % VTK_TETRA
      lf   = [ 1, 2, 4 ;   % face opposite node 3
               2, 3, 4 ;   % face opposite node 1
               3, 1, 4 ;   % face opposite node 2
               1, 3, 2 ];  % base face (opposite node 4)
      nNPF = 3;

    case {12, 29}    % VTK_HEXAHEDRON  /  VTK_TRIQUADRATIC_HEXAHEDRON
      lf   = [ 1, 4, 3, 2 ;   % bottom  (z-)
               5, 6, 7, 8 ;   % top     (z+)
               1, 2, 6, 5 ;   % front   (y-)
               2, 3, 7, 6 ;   % right   (x+)
               3, 4, 8, 7 ;   % back    (y+)
               4, 1, 5, 8 ];  % left    (x-)
      nNPF = 4;

    otherwise
      error('processFaces:unsupportedVTKType', ...
           ['Unsupported VTK 3-D cell type: %d.\n'     ...
            'Supported types: 10 (tet4), 12 (hex8), ' ...
            '29 (hex27).'], vtkId);
  end
end


% =========================================================================
%  numCorners  —  number of corner nodes for each supported VTK cell type
% =========================================================================
function n = numCorners(vtkId)
  switch vtkId
    case 10,       n = 4;   % tet4:  all 4 nodes are corners
    case {12, 29}, n = 8;   % hex8 / hex27:  first 8 nodes are corners
    otherwise
      error('processFaces:unsupportedVTKType', ...
            'numCorners: unsupported VTK type %d.', vtkId);
  end
end

