function processGeometry(grid)
% process alla geometric informations of the grid starting from minimal
% connectivity informations

% cells
cells = grid.cells;
surf = grid.surfaces;

cells.num = size(cells.connectivity,1);
surf.num = size(surf.connectivity,1); 

if cells.num > 0
  grid.nDim = 3;
else
  grid.nDim = 2;
end

% check if mixed cell shapes are present
vtkTypes = reshape(unique(cells.VTKType),1,[]);
if numel(vtkTypes) > 1
  grid.isMixed = true;
end

nc = cells.num;
ns = surfaces.num;

cells.volume = zeros(nc,1);
cells.center = zeros(nc,3);

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



function processShape(grid,vtk3d,cells,surf)

i = find(grid.vtkType(:,2)==vtk3d);
vtk2d =  grid.vtkType(i,1);

volShape = FiniteElementType.create(vtk3d,grid);
[vols, cellCenters] = getSizeAndCentroid(volShape);

% process grid.faces.connectivity and grid.faces.neighbors
processFaces(grid,vtk3d);

surfShape = FiniteElementType.create(vtk2d,grid);
[areas, faceCenters] = getSizeAndCentroid(surfShape);


end





