function processGeometry(grid)
% process alla geometric informations of the grid starting from minimal
% connectivity informations

if grid.isProcessed
  warning('The grid geometry has already been processed, geometry will be restored')
end

initializeGrid(grid);

for vtkId = grid.cells.vtkTypes
  processShape(grid,vtkId);
end


% fix normal orientation and finalize surface geometry
fixNormals(grid);

fId = grid.surfaces.faceId;
grid.surfaces.area = grid.faces.area(fId);
grid.surfaces.center = grid.faces.center(fId,:);
grid.surfaces.normal = grid.faces.normal(fId,:);


% set grid processed
grid.isProcessed = true;

end



function processShape(grid,vtk3d)

idC = grid.getCellsByVTKId(vtk3d);

% process faces
grid.processFaces(vtk3d);

% compute geometry of the cells
volShape = FiniteElementType.create(vtk3d,grid);
[grid.cells.volume(idC), grid.cells.center(idC,:)] = getSizeAndCentroid(volShape);

end


function fixNormals(grid)

% fix the orientation of normals in the grid
neigh = grid.faces.neighbors(:,1);
d = grid.faces.center - grid.cells.center(neigh,:);

toFix = sum(d.*grid.faces.normal,2) < 0;
grid.faces.normal(toFix,:) = - grid.faces.normal(toFix,:);



end






