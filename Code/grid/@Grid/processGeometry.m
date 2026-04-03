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
  [cells,faces] = processShape(vtkId,cells);
end

% link faces to surfaces
surf.faceId = []; % this can be used by getSurface to derive all properties of related face


grid.cells = cells;
grid.surfaces = surf;
grid.faces = faces;

end



function [cells,faces] = processShape(vtkId,cells)

switch vtkId
  case 9
    [cells,faces] = processTetra(cells);
  case 12
    [cells,faces] = processHexa(cells);
  case 29
    gresLog().warning(2,"Faces are not supported with quadratic hexahedra. Finite volume based discretization are cannot be used.")
    [cells,faces] = processHexaQuad(cells);
  otherwise
    error('VTK type %i is not yet supported');
end

end



function processTetra(obj)

% use tetra utils
tetra = Tetrahedron(obj);

% make topology
[obj.cells.volume, obj.cells.center] = findVolumeAndCentroid(tetra);



end


function processHexa(cells)

% use hexa utils
hexa = Hexahedron();

end

