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
      [~,geometricNodes] = ismember(matrix.coordinates, ...
        unique(matrix.coordinates,'rows'),'rows');
      for k = 1:2
        fracture = fractures{k};
        testCase.verifyEqual([fracture.nDim fracture.nCoord],[2 3]);
        testCase.verifyGreaterThan(fracture.surfaces.num,0);
        pairs = fracture.cells.matrixFaces;
        testCase.verifyGreaterThan(fracture.cells.num,0);
        fractureCells = fracture.cells.loc2glob(fracture.cells.connectivity);
        expected = sort(geometricNodes(fractureCells),2);
        testCase.verifyEqual(sort(geometricNodes(matrixFaces(pairs(:,1),:)),2),expected);
        testCase.verifyEqual(sort(geometricNodes(matrixFaces(pairs(:,2),:)),2),expected);
        testCase.verifyEqual(matrix.faces.neighbors(pairs(:),2),zeros(numel(pairs),1));
        owners = matrix.faces.neighbors(pairs(:),1);
        testCase.verifyTrue(all(any(matrixCellFaces(owners,:) == pairs(:),2)));
        testCase.verifyDisjointOwners(matrix,pairs);
        testCase.verifyEqual(matrix.faces.normal(pairs(:,1),:), ...
          -matrix.faces.normal(pairs(:,2),:));
        testCase.verifyEqual(sum(fracture.cells.area),.36,'AbsTol',1e-10);
        fractureEdges{k} = sort([fractureCells(:,[1 2]); fractureCells(:,[2 3]); ...
          fractureCells(:,[3 1])],2);
      end
      testCase.verifyNotEmpty(intersect(fractureEdges{1},fractureEdges{2},'rows'));
      testCase.verifyEmpty(intersect(matrix.surfaces.faceId, ...
        [fractures{1}.cells.matrixFaces(:); fractures{2}.cells.matrixFaces(:)]));
    end

    function buildsStructuredFractureGrid(testCase)
      bounds = [0 1; 0 1; 0 1];
      perimeter = [.5 .2 .2; .5 .8 .2; .5 .8 .8; .5 .2 .8];
      [matrix,fractures] = createFracturedGrid(bounds,perimeter,struct('structured',true));
      fracture = fractures{1};
      testCase.verifyTrue(matrix.isStructured());
      testCase.verifyTrue(all(matrix.cells.VTKType == VTKType.Hexa));
      testCase.verifyTrue(all(fracture.cells.VTKType == VTKType.Quad));
      testCase.verifyEqual(sum(matrix.cells.volume),1,'AbsTol',1e-12);
      testCase.verifyEqual(sum(fracture.cells.area),.36,'AbsTol',1e-12);
      testCase.verifyDisjointOwners(matrix,fracture.cells.matrixFaces);

      perimeter(:,1) = perimeter(:,2);
      testCase.verifyError(@() createFracturedGrid(bounds,perimeter,'structured',true), ...
        'GReS:StructuredFractureAlignment');
    end
  end

  methods (Access = private)
    function verifyDisjointOwners(testCase,matrix,pairs)
      owners = reshape(matrix.faces.neighbors(pairs(:),1),size(pairs));
      left = matrix.cells.connectivity(owners(:,1),:);
      right = matrix.cells.connectivity(owners(:,2),:);
      testCase.verifyFalse(any(left == permute(right,[1 3 2]),'all'));
    end
  end
end
