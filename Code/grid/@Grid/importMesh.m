function importMesh(grid, fileName)
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
    [grid.coordinates, elems] = mxImportVTKmesh(char(fileName));
    elems = double(elems);

  case 'msh'
    [grid.coordinates, elems, ~] = mxImportGMSHmesh(char(fileName));
    % the third input contains region names and is currently not supported
    % by GReS
    elems = double(elems);
    elems = gmshToVTK(elems);
  otherwise
    error("Unsupported mesh format '.%s'. Valid formats are: 'vtk', 'msh'.", extension);
end

% -------------------------------------------------
% DIMENSIONALITY
% -------------------------------------------------
if any(grid.coordinates(:,3) ~= 0)
  grid.nDim = 3;
else
  grid.nDim = 2;
end
grid.nNodes = size(grid.coordinates, 1);

% -------------------------------------------------
% 3D CELL DATA
% -------------------------------------------------
if any(~ismember(elems(:,1), grid.vtkType))
  error(['Unsupported  elements in the mesh.\n', ...
    'Supported types: 3-node triangles (VTKType = 5), ', ...
    '4-node quadrilaterals (VTKType = 9).',...
    '9-node quadrilaterals (VTKType = 28).',...
    'Supported 3D types: 4-node tetrahedra (VTKType = 10), ', ...
    '8-node hexahedra (VTKType = 12).',...
    '27-node hexahedra (VTKType = 29).']);
end

ID = ismember(elems(:,1), grid.vtkType(:,2));
grid.cells.numVerts  = elems(ID, 3);
nVerts            = max(grid.cells.numVerts);
grid.cells.connectivity = Grid.makeConnectivity(elems(ID, 4:nVerts+3));
grid.cells.VTKType   = elems(ID, 1);
grid.cells.tag       = elems(ID, 2);

if all(grid.cells.tag == 0)
  grid.cells.tag  = grid.cells.tag + 1;
end

% -------------------------------------------------
% 2D SURFACE DATA - tagged boundary surfaces (for bc and grid coupling)
% -------------------------------------------------

ID = ismember(elems(:,1), grid.vtkType(:,1));
grid.surfaces.numVerts = elems(ID, 3);
nVerts = max(grid.surfaces.numVerts);
grid.surfaces.connectivity = Grid.makeConnectivity(elems(ID, 4:nVerts+3));
grid.surfaces.VTKType  = elems(ID, 1);
grid.surfaces.tag      = elems(ID, 2);

% -------------------------------------------------
% 1D EDGE DATA
% -------------------------------------------------

% for now, edge data is derived from external surfaces only
% potentially, this is also a field of the imported grid

% ID = ismember(elems(:,1), grid.edgeVTK);
% grid.edgeNumVerts = elems(ID, 3);
% nVerts           = max(grid.edgeNumVerts);
% grid.edges        = elems(ID, 4:nVerts+3);
% grid.edgeTag      = elems(ID, 2);
% grid.edgeVTKType  = elems(ID, 1);
% grid.nEdges       = length(grid.edgeTag);

processGeometry(grid);

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
%     GMSH 10 → VTK  28  (9-node quad → stored as 4-node quad)
%     GMSH 12 → VTK 29  (27-node hex → stored as 8-node hex)

gmshIDs = [2  3,  4,  10,  5,  12];
vtkIDs  = [5,  9, 10,   28, 12,  29];

newCol = elems(:,1);
for k = 1:length(gmshIDs)
  newCol(elems(:,1) == gmshIDs(k)) = vtkIDs(k);
end
elems(:,1) = newCol;

end
