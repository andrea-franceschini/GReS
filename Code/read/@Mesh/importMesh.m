function importMesh(Mesh, fileName)
% IMPORTMESH  Import a mesh file and populate the GReS Mesh object
% properties.
%
%   IMPORTMESH(mesh, fileName) reads mesh data from the specified file and
%   stores element connectivity, node coordinates, tags, and region data
%   in the Mesh object mesh.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS
% -------------------------------------------------------------------------
%
%   fileName  - (string) Path to the mesh file. The file extension
%               determines the format:
%                 '.vtk'  — VTK unstructured grid (via mxImportVTKmesh)
%                 '.msh'  — GMSH mesh (via mxImportGMSHmesh)
%
% -------------------------------------------------------------------------
% SUPPORTED ELEMENT TYPES
% -------------------------------------------------------------------------
%
%   3D cells:    4-node tetrahedra  (VTKType = 10)
%                8-node hexahedra   (VTKType = 12)
%                27-node hexahedra  (VTKType = 29)
%
%   2D surfaces: 3-node triangles        (VTKType = 5)
%                4-node quadrilaterals   (VTKType = 9)
%                9-node quadrilaterlas   (VTKTYPE = 28)
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%
%   - If no physical volume tag is defined in the mesh, all cells are
%     assigned tag 1 by default.
%   - Region data (named physical groups) is populated only when available.
%     GMSH files always provide regions; VTK files provide them only if
%     the 'regions' variable is present in the imported data.
%
% -------------------------------------------------------------------------
% SEE ALSO
% -------------------------------------------------------------------------
%   mxImportVTKmesh, mxImportGMSHmesh, gmshToVTK

str = split(fileName, '.');
extension = str{2};

switch extension
  case 'vtk'
    [Mesh.coordinates, elems] = mxImportVTKmesh(char(fileName));
    elems = double(elems);
    regions = [];
    Mesh.meshType = 'Unstructured';

  case 'msh'
    [Mesh.coordinates, elems, ~] = mxImportGMSHmesh(char(fileName));
    % the third input contains region names and is currently not supported
    % by GReS
    elems = double(elems);
    elems = gmshToVTK(elems);
    Mesh.meshType = 'Unstructured';

  otherwise
    error("Unsupported mesh format '.%s'. Valid formats are: 'vtk', 'msh'.", extension);
end

% -------------------------------------------------
% DIMENSIONALITY
% -------------------------------------------------
if any(Mesh.coordinates(:,3) ~= 0)
  Mesh.nDim = 3;
else
  Mesh.nDim = 2;
end
Mesh.nNodes = size(Mesh.coordinates, 1);

% -------------------------------------------------
% 3D CELL DATA
% -------------------------------------------------
ID = ismember(elems(:,1), Mesh.cellVTK);
Mesh.cellNumVerts  = elems(ID, 3);
nVerts            = max(Mesh.cellNumVerts);
Mesh.cells         = elems(ID, 4:nVerts+3);
Mesh.cellVTKType   = elems(ID, 1);
Mesh.cellTag       = elems(ID, 2);
Mesh.nCells        = length(Mesh.cellTag);

if any(~ismember(Mesh.cellVTKType, [10, 12, 29]))
  error(['Unsupported 3D elements in the mesh.\n', ...
    'Supported types: 4-node tetrahedra (VTKType = 10), ', ...
    '8-node hexahedra (VTKType = 12).',...
    '27-node hexahedra (VTKType = 29).']);
end

if all(Mesh.cellTag == 0)
  Mesh.cellTag  = Mesh.cellTag + 1;
  Mesh.nCellTag = 1;
else
  Mesh.nCellTag = max(Mesh.cellTag);
end

% -------------------------------------------------
% 2D SURFACE DATA
% -------------------------------------------------
ID = ismember(elems(:,1), Mesh.surfaceVTK);
Mesh.surfaceNumVerts = elems(ID, 3);
nVerts              = max(Mesh.surfaceNumVerts);
Mesh.surfaces        = elems(ID, 4:nVerts+3);
Mesh.surfaceVTKType  = elems(ID, 1);
Mesh.surfaceTag      = elems(ID, 2);
Mesh.nSurfaces       = length(Mesh.surfaceTag);
Mesh.nSurfaceTag     = max(Mesh.surfaceTag);

if any(~ismember(Mesh.surfaceVTKType, [5, 9, 28]))
  error(['Unsupported 2D elements in the mesh.\n', ...
    'Supported types: 3-node triangles (VTKType = 5), ', ...
    '4-node quadrilaterals (VTKType = 9).',...
    '9-node quadrilaterals (VTKType = 28).']);
end

% -------------------------------------------------
% 1D EDGE DATA
% -------------------------------------------------
ID = ismember(elems(:,1), Mesh.edgeVTK);
Mesh.edgeNumVerts = elems(ID, 3);
nVerts           = max(Mesh.edgeNumVerts);
Mesh.edges        = elems(ID, 4:nVerts+3);
Mesh.edgeTag      = elems(ID, 2);
Mesh.edgeVTKType  = elems(ID, 1);
Mesh.nEdges       = length(Mesh.edgeTag);

end


% =========================================================================

function elems = gmshToVTK(elems)
% GMSHTOVTK  Convert GMSH element table to VTK-compatible format.
%
%   elems = GMSHTOVTK(obj, elems) remaps GMSH element type IDs to their
%   VTK equivalents so that the unified importMesh storage block can
%   process both VTK and GMSH meshes identically.
%
%   The input 'elems' matrix uses GMSH type IDs in column 1.
%   The output 'elems' matrix uses VTK type IDs in column 1.
%
%   GMSH-to-VTK type mapping:
%     GMSH  2 → VTK  5  (3-node triangle)
%     GMSH  3 → VTK  9  (4-node quadrilateral)
%     GMSH  4 → VTK 10  (4-node tetrahedron)
%     GMSH  5 → VTK 12  (8-node hexahedron)
%     GMSH 10 → VTK  9  (9-node quad → stored as 4-node quad)
%     GMSH 12 → VTK 12  (27-node hex → stored as 8-node hex)

gmshIDs = [2,  3,  4,  10,  5,  12];
vtkIDs  = [5,  9, 10,   28, 12,  29];

newCol = elems(:,1);
for k = 1:length(gmshIDs)
  newCol(elems(:,1) == gmshIDs(k)) = vtkIDs(k);
end
elems(:,1) = newCol;

end
