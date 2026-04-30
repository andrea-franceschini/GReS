function initializeGrid(grid)

grid.nNodes = size(grid.coordinates,1);
grid.cells.num = size(grid.cells.connectivity,1);
grid.surfaces.num = size(grid.surfaces.connectivity,1);
nc = grid.cells.num;
ns = grid.surfaces.num;

grid.cells.cells2faces = ArrayOfArrays();
grid.cells.cells2localFaces = ArrayOfArrays();
%grid.surfaces.faceId = [];


if grid.nNodes > 0    % grid already populated


  if nc > 0
    grid.cells.nTag = numel(unique(grid.cells.tag));
    % check if mixed cell shapes are present
    vtkTypes = reshape(unique(grid.cells.VTKType),1,[]);
    grid.cells.vtkTypes = vtkTypes;
    
    if numel(vtkTypes) > 1
      grid.isMixed = true;
    end

  end

  if ns > 0
    grid.surfaces.vtkTypes = reshape(unique(grid.surfaces.VTKType),1,[]);
    grid.surfaces.nTag = numel(unique(grid.surfaces.tag));
  end

  if grid.cells.num > 0
    grid.nDim = 3;
  elseif grid.surfaces.num > 0
    grid.nDim = 2;
  end

  return

end

% grid must be populated

grid.cells.volume = zeros(nc,1);
grid.cells.center = zeros(nc,3);
grid.cells.VTKType = zeros(nc,1);
grid.cells.tag = zeros(nc,1);
grid.cells.nTag = 0;
grid.cells.numVerts = zeros(nc,1);


grid.surfaces.area = zeros(ns,1);
grid.surfaces.numVerts = zeros(ns,1);
grid.surfaces.VTKType = zeros(ns,1);
grid.surfaces.tag = zeros(ns,1);
grid.surfaces.nTag = 0;
grid.surfaces.center = zeros(ns,3);
grid.surfaces.normal = zeros(ns,3);
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

