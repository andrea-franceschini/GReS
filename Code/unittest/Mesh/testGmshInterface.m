classdef testGmshInterface < matlab.unittest.TestCase
  
  methods (Test)
    function splitsIntersectingFractures(testCase)
      testCase.assumeTrue(gresgmsh.isAvailable(), ...
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
        pairs = fracture.cells.matrixFaces;
        testCase.verifyGreaterThan(fracture.cells.num,0);
        fractureCells = fracture.cells.loc2glob(fracture.cells.connectivity);
        expected = sort(geometricNodes(fractureCells),2);
        testCase.verifyEqual(sort(geometricNodes(matrixFaces(pairs(:,1),:)),2),expected);
        testCase.verifyEqual(sort(geometricNodes(matrixFaces(pairs(:,2),:)),2),expected);
        testCase.verifyEqual(matrix.faces.neighbors(pairs(:),2),zeros(numel(pairs),1));
        owners = matrix.faces.neighbors(pairs(:),1);
        testCase.verifyTrue(all(any(matrixCellFaces(owners,:) == pairs(:),2)));
        testCase.verifyIntactFaceNodes(matrix);
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
      testCase.verifyClass(matrix,'StructuredGrid');
      testCase.verifyTrue(matrix.isStructured());
      testCase.verifyTrue(all(matrix.cells.VTKType == VTKType.Hexa));
      testCase.verifyTrue(all(fracture.cells.VTKType == VTKType.Quad));
      testCase.verifyEqual(sum(matrix.cells.volume),1,'AbsTol',1e-12);
      testCase.verifyEqual(sum(fracture.cells.area),.36,'AbsTol',1e-12);
      testCase.verifyIntactFaceNodes(matrix);

      perimeter(:,1) = perimeter(:,2);
      testCase.verifyError(@() createFracturedGrid(bounds,perimeter,'structured',true), ...
        'GReS:StructuredFractureAlignment');
    end

    function buildsIntersectingStructuredFractures(testCase)
      bounds = [0 1; 0 1; 0 1];
      fractures = {struct('center',[.5 .5 .5],'normal',[1 0 0],'extents',[.6 .6]), ...
        struct('center',[.5 .5 .5],'normal',[0 -1 0],'extents',[.6 .6]), ...
        struct('center',[.5 .5 .5],'normal',[0 0 1],'extents',[.6 .6])};
      [matrix,grids] = createFracturedGrid(bounds,fractures,'structured',true, ...
        'matrixMeshSize',.25,'fractureMeshSize',.15);
      testCase.verifyClass(matrix,'StructuredGrid');
      atCrossing = all(abs(matrix.coordinates - [.5 .5 .5]) < 1e-12,2);
      testCase.verifyEqual(nnz(atCrossing),8);
      testCase.verifyEqual(sum(matrix.cells.volume),1,'AbsTol',1e-12);
      for k = 1:numel(grids)
        grid = grids{k};
        testCase.verifyEqual(sum(grid.cells.area),.36,'AbsTol',1e-12);
        testCase.verifyIntactFaceNodes(matrix);
        testCase.verifyGreaterThan(grid.cells.normal * fractures{k}.normal.', ...
          zeros(grid.cells.num,1));
      end
    end

    function keepsStructuredTipNodesShared(testCase)
      testCase.verifyTipNodes(true);
    end

    function keepsGmshTipNodesShared(testCase)
      testCase.assumeTrue(gresgmsh.isAvailable(), ...
        'Optional Gmsh MEX is unavailable.');
      testCase.verifyTipNodes(false);
    end

    function buildsStructuredMeshes(testCase)
      grids = {structuredMesh([0 .2 1],[0 .5 1]), ...
        structuredMesh(2,2,[0 1],[0 1],'triangle',true), ...
        structuredMesh(2,2,2,[0 1],[0 1],[0 1]), ...
        structuredMesh(2,2,2,[0 1],[0 1],[0 1],'tetra',true)};
      for k = 1:numel(grids)
        testCase.verifyClass(grids{k},'StructuredGrid');
        grids{k}.processGeometry();
        testCase.verifyEqual(sum(grids{k}.cells.volume),1,'AbsTol',1e-12);
      end
      [matrix,fractures] = createFracturedGrid([0 1;0 1;0 1],[], ...
        'structured',true,'meshSize',.5);
      testCase.verifyClass(matrix,'StructuredGrid');
      testCase.verifyEmpty(fractures);
      testCase.verifyEqual(matrix.cells.num,8);
    end
  end

  methods (Access = private)
    function verifyIntactFaceNodes(testCase,matrix)
      internal = ~matrix.faces.isBoundary;
      owners = matrix.faces.neighbors(internal,:);
      faceNodes = matrix.faces.connectivity.toMatrix();
      faceNodes = faceNodes(internal,:);
      for side = 1:2
        cellNodes = matrix.cells.connectivity(owners(:,side),:);
        matched = any(faceNodes == permute(cellNodes,[1 3 2]),3);
        testCase.verifyTrue(all(matched(faceNodes > 0)));
      end
    end

    function verifyTipNodes(testCase,structured)
      bounds = [0 1;0 1;0 1];
      polygon = [.5 .2 .2;.5 .8 .2;.5 .8 .8;.5 .2 .8];
      [matrix,fractures] = createFracturedGrid(bounds,polygon, ...
        'structured',structured,'meshSize',.1);
      fracture = fractures{1};
      boundaryEdges = fracture.faces.connectivity.toMatrix();
      tips = unique(boundaryEdges(fracture.faces.isBoundary,:));
      tips(tips == 0) = [];
      [xyz,~,node] = unique(matrix.coordinates,'rows');
      copies = accumarray(node,1);
      [found,geometricNode] = ismember(fracture.coordinates,xyz,'rows');
      testCase.verifyTrue(all(found));
      testCase.verifyEqual(copies(geometricNode(tips)),ones(numel(tips),1));
      interior = setdiff((1:fracture.nNodes).',tips);
      testCase.verifyNotEmpty(interior);
      testCase.verifyEqual(copies(geometricNode(interior)),2*ones(numel(interior),1));
      testCase.verifyIntactFaceNodes(matrix);
    end
  end
end
