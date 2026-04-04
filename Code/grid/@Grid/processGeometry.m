function processGeometry(grid)
% process alla geometric informations of the grid starting from minimal
% connectivity informations

% cells
cells = grid.cells;
surf = grid.surfaces;

cells.num = size(cells.connectivity,1);
surf.num = size(surf.connectivity,1); 

% check if mixed cell shapes are present
vtkTypes = reshape(unique(cells.VTKtype),1,[]);
if numel(vtkTypes) > 1
  grid.isMixed = true;
end

% process vtk types
for vtkId = vtkTypes
  % return face connectivity and 3D-2D geometry for given vtk type
  [cells,faces,surf] = processShape(grid,vtkId,cells,surf);
end

% link faces to surfaces
surf.faceId = []; % this can be used by getSurface to derive all properties of related face


grid.cells = cells;
grid.surfaces = surf;
grid.faces = faces;

end



function [cells,faces] = processShape(grid,vtk3d,cells,surf)

i = find(grid.vtkType(:,2)==vtk3d);
vtk2d =  obj.vtkId(i,1);

volShape = FiniteElementType.create(vtk3d,grid);
[vols, cellCenters] = getSizeAndCentroid(volShape);

% process grid.faces.connectivity and grid.faces.neighbors
processFaces(obj,vtk3d);

surfShape = FiniteElementType.create(vtk2d,grid);
[areas, faceCenters] = getSizeAndCentroid(surfShape);


end





