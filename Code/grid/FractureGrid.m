classdef FractureGrid < Grid
  % A triangular grid coupled to the two matrix faces at each triangle.

  properties
    matrix
  end

  methods
    function obj = FractureGrid(matrix, triangles, matrixFaces, fractureId)
      obj.matrix = matrix;
      normal = matrix.faces.normal(matrixFaces(:,1),:);
      xyz = matrix.coordinates;
      winding = cross(xyz(triangles(:,2),:) - xyz(triangles(:,1),:), ...
                      xyz(triangles(:,3),:) - xyz(triangles(:,1),:),2);
      flip = dot(winding,normal,2) < 0;
      triangles(flip,[2 3]) = triangles(flip,[3 2]);
      [nodes,~,localNodes] = unique(triangles(:),'stable');
      obj.coordinates = matrix.coordinates(nodes,:);
      obj.nNodes = numel(nodes);
      obj.nDim = 2;
      obj.cells.connectivity = reshape(localNodes,size(triangles));
      obj.cells.num = size(triangles,1);
      obj.cells.numVerts = repmat(3,obj.cells.num,1);
      obj.cells.tag = repmat(fractureId,obj.cells.num,1);
      obj.cells.nTag = 1;
      obj.cells.vtkTypes = double(VTKType.Tri);
      obj.cells.VTKType = repmat(obj.cells.vtkTypes,obj.cells.num,1);
      obj.cells.loc2glob = nodes;
      obj.cells.matrixFaces = matrixFaces;
      obj.cells.center = matrix.faces.center(matrixFaces(:,1),:);
      obj.cells.area = matrix.faces.area(matrixFaces(:,1));
      obj.cells.normal = normal;
    end
  end
end
