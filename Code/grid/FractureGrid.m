classdef FractureGrid < Grid
  % A surface grid coupled to two matrix faces per cell.

  properties
    matrix
  end

  methods
    function obj = FractureGrid(matrix, connectivity, matrixFaces, fractureId)


      obj.matrix = matrix;
      normal = matrix.faces.normal(matrixFaces(:,1),:);
      xyz = matrix.coordinates;
      nVerts = size(connectivity,2);
      winding = cross(xyz(connectivity(:,2),:) - xyz(connectivity(:,1),:), ...
        xyz(connectivity(:,3),:) - xyz(connectivity(:,1),:),2);
      flip = dot(winding,normal,2) < 0;
      connectivity(flip,2:end) = connectivity(flip,end:-1:2);
      [nodes,~,localNodes] = unique(connectivity(:),'stable');
      obj.coordinates = matrix.coordinates(nodes,:);
      obj.nNodes = numel(nodes);
      obj.nDim = 2;
      obj.cells.connectivity = reshape(localNodes,size(connectivity));
      obj.cells.num = size(connectivity,1);
      obj.cells.numVerts = repmat(nVerts,obj.cells.num,1);
      obj.cells.tag = repmat(fractureId,obj.cells.num,1);
      obj.cells.nTag = 1;
      obj.cells.vtkTypes = double(matrix.faces.VTKType(matrixFaces(1,1)));
      obj.cells.VTKType = repmat(obj.cells.vtkTypes,obj.cells.num,1);
      obj.cells.loc2glob = nodes;
      obj.cells.matrixFaces = matrixFaces;
      obj.cells.center = matrix.faces.center(matrixFaces(:,1),:);
      obj.cells.area = matrix.faces.area(matrixFaces(:,1));
      obj.cells.normal = normal;
      obj.processGeometry();
    end
  end
end
