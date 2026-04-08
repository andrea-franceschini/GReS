function processGeometry(grid)
% process alla geometric informations of the grid starting from minimal
% connectivity informations

initializeGrid(grid);

% process vtk types
grid.cells.vtkTypes = reshape(unique(grid.cells.VTKType),1,[]);
grid.surfaces.vtkTypes = reshape(unique(grid.surfaces.VTKType),1,[]);

for vtkId = grid.cells.vtkTypes
  processShape(grid,vtkId);
end


% fix normal orientation and finalize surface geometry
fixNormals(grid);

fId = grid.surfaces.faceId;
grid.surfaces.area = grid.faces.area(fId);
grid.surfaces.center = grid.faces.center(fId,:);
grid.surfaces.normal = grid.faces.normal(fId,:);

end



function processShape(grid,vtk3d)

idC = grid.getCellsByVTKId(vtk3d);

% process faces
grid.processFaces(vtk3d);

% compute geometry of the cells
volShape = FiniteElementType.create(vtk3d,grid);
[grid.cells.volume(idC), grid.cells.center(idC,:)] = getSizeAndCentroid(volShape);

end



function initializeGrid(grid)

grid.nNodes = size(grid.coordinates,1);
grid.cells.num = size(grid.cells.connectivity,1);
grid.surfaces.num = size(grid.surfaces.connectivity,1); 
nc = grid.cells.num;
ns = grid.surfaces.num;

if grid.cells.num > 0
  grid.nDim = 3;
else
  grid.nDim = 2;
end

% check if mixed cell shapes are present
vtkTypes = reshape(unique(grid.cells.VTKType),1,[]);
if numel(vtkTypes) > 1
  grid.isMixed = true;
end


grid.cells.volume = zeros(nc,1);
grid.cells.center = zeros(nc,3);
grid.cells.cells2faces = ArrayOfArrays();
grid.cells.cells2localFaces = ArrayOfArrays();



grid.surfaces.area = zeros(ns,1);
grid.surfaces.center = zeros(ns,3);
grid.surfaces.faceId = zeros(ns,1);

grid.faces.num = 0;
grid.faces.numVerts = zeros(0,1);
grid.faces.area = zeros(0,1);
grid.faces.normal = zeros(0,3);
grid.faces.center = zeros(0,3);
grid.faces.neighbors = zeros(0,2);
grid.faces.isBoundary = false(0,1);
grid.faces.connectivity = ArrayOfArrays();

end



function fixNormals(grid)

% fix the orientation of normals in the grid
neigh = grid.faces.neighbors(:,1);
d = grid.faces.center - grid.cells.center(neigh,:);

toFix = sum(d.*grid.faces.normal,2) < 0;
grid.faces.normal(toFix,:) = - grid.faces.normal(toFix,:);



end






