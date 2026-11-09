classdef testGmshInterface < matlab.unittest.TestCase
  methods (Test)
    function splitsIntersectingFractures(testCase)
      testCase.assumeEqual(exist('gmshMex','file'),3, ...
        'Optional Gmsh MEX is unavailable.');
      bounds = [0 1; 0 1; 0 1];
      perimeter = [.5 .2 .2; .5 .8 .2; .5 .8 .8; .5 .2 .8];
      rectangle = struct('center',[.5 .5 .5],'normal',[0 1 0],'extents',[.6 .6]);
      [matrix, fractures] = createFracturedGrid(bounds,{perimeter,rectangle}, ...
        'matrixMeshSize',.24,'fractureMeshSize',.10,'transitionFactor',1.5);
      testCase.verifyNumElements(fractures,2);
      matrixFaces = matrix.faces.connectivity.toMatrix();
      matrixCellFaces = matrix.cells.cells2faces.toMatrix();
      fractureEdges = cell(1,2);
      for k = 1:2
        fracture = fractures{k};
        testCase.verifyEqual([fracture.nDim fracture.nCoord],[2 3]);
        testCase.verifyGreaterThan(fracture.surfaces.num,0);
        pairs = fracture.cells.matrixFaces;
        testCase.verifyGreaterThan(fracture.cells.num,0);
        triangles = fracture.cells.loc2glob(fracture.cells.connectivity);
        testCase.verifyEqual(matrixFaces(pairs(:,1),:),sort(triangles,2));
        testCase.verifyEqual(matrixFaces(pairs(:,1),:),matrixFaces(pairs(:,2),:));
        testCase.verifyEqual(matrix.faces.neighbors(pairs(:),2),zeros(numel(pairs),1));
        owners = matrix.faces.neighbors(pairs(:),1);
        testCase.verifyTrue(all(any(matrixCellFaces(owners,:) == pairs(:),2)));
        testCase.verifyEqual(matrix.faces.normal(pairs(:,1),:), ...
          -matrix.faces.normal(pairs(:,2),:));
        testCase.verifyEqual(sum(fracture.cells.area),.36,'AbsTol',1e-10);
        fractureEdges{k} = sort([triangles(:,[1 2]); triangles(:,[2 3]); ...
          triangles(:,[3 1])],2);
      end
      testCase.verifyNotEmpty(intersect(fractureEdges{1},fractureEdges{2},'rows'));
      testCase.verifyEmpty(intersect(matrix.surfaces.faceId, ...
        [fractures{1}.cells.matrixFaces(:); fractures{2}.cells.matrixFaces(:)]));
    end
  end
end
