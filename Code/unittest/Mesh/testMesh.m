classdef testMesh < matlab.unittest.TestCase

   methods(Test)
    % Test methods
    function testReadVTK(testCase)
      grid = Grid();
      grid.importMesh(fullfile(fileparts(mfilename('fullpath')),'testMesh.vtk'));
      grid.processGeometry;
      testCase.verifyEqual(grid.nNodes,919);
      testCase.verifyEqual(grid.nDim,3);
      testCase.verifyEqual(grid.cells.num,560);
      testCase.verifyEqual(grid.surfaces.num,624);
    end

    function testReadQuadratic2DMesh(testCase)
      grid = Grid();
      grid.importMesh(fullfile(fileparts(mfilename('fullpath')),'quad9.msh'));
      grid.processGeometry();
      testCase.verifyEqual([grid.nDim grid.nCoord],[2 3]);
      testCase.verifyEqual(grid.cells.num,1);
      testCase.verifyEqual(grid.getCellNodes(1),1:9);
      testCase.verifyEqual(grid.cells.volume,1,'AbsTol',1e-12);
      testCase.verifyEqual(grid.surfaces.tag,(1:4).');
    end

    function testPlanarGridAndBoundaryExtraction(testCase)
      grid = structuredMesh(2,2,[0 1],[0 1]);
      grid.processGeometry();
      testCase.verifyEqual([grid.nDim grid.nCoord],[2 2]);
      testCase.verifyEqual(grid.cells.volume,repmat(.25,4,1),'AbsTol',1e-12);
      testCase.verifyEqual(grid.faces.num,12);
      testCase.verifyEqual(grid.surfaces.num,8);
      testCase.verifyEqual(unique(grid.surfaces.VTKType),double(VTKType.Line));
      boundary = grid.getSurfaceGrid(1);
      testCase.verifyEqual(boundary.nDim,1);
      testCase.verifyEqual(boundary.cells.num,2);
      testCase.verifyEqual(boundary.cells.faceId,grid.surfaces.faceId(grid.surfaces.tag == 1));
      [~,influence] = grid.getNodeInfluence(entityField.surface,1);
      testCase.verifyEqual(influence,[.25;.25],'AbsTol',1e-12);
      grid.processGeometry();
      testCase.verifyEqual([grid.nDim grid.nCoord],[2 2]);
    end

    function testMixedEmbeddedCells(testCase)
      grid = Grid();
      grid.coordinates = [0 0 0; 1 0 1; 1 1 1; 0 1 0; 2 0 2];
      grid.cells.connectivity = [1 2 3 4; 2 5 3 0];
      grid.cells.VTKType = [9;5];
      grid.cells.tag = [1;1];
      grid.cells.numVerts = [4;3];
      grid.processGeometry();
      testCase.verifyEqual([grid.nDim grid.nCoord],[2 3]);
      testCase.verifyEqual(grid.cells.volume,sqrt(2)*[1;.5],'AbsTol',1e-12);
      testCase.verifyEqual(grid.faces.num,6);
      testCase.verifyEqual(grid.surfaces.num,5);
      internal = find(~grid.faces.isBoundary);
      testCase.verifyNumElements(internal,1);
      for k = 1:2
        testCase.verifyTrue(ismember(internal,grid.cells.cells2faces.getRowsMatrix(k)));
      end
    end

    function testExtractedSurfaceCells(testCase)
      parent = structuredMesh(2,2,1,[0 1],[0 1],[0 1]);
      parent.processGeometry();
      patch = parent.getSurfaceGrid(2);
      testCase.verifyEqual([patch.nDim patch.nCoord],[2 3]);
      testCase.verifyEqual(patch.cells.num,4);
      testCase.verifyEqual(patch.surfaces.num,8);
      testCase.verifyEqual(sum(patch.cells.volume),1,'AbsTol',1e-12);
      subset = patch.getCellGrid([true;false;false;false]);
      subset.processGeometry();
      testCase.verifyEqual(subset.cells.faceId,patch.cells.faceId(1));
      testCase.verifyEqual(subset.cells.normal,parent.faces.normal(subset.cells.faceId,:));
      nodes = subset.cells.loc2glob(subset.getCellNodes(1));
      parentNodes = parent.faces.connectivity.getRowsMatrix(subset.cells.faceId);
      testCase.verifyEqual(sort(nodes(:)),sort(parentNodes(:)));
    end

    function testStructuredGrids(testCase)
      grid = structuredMesh(5,5,5,[0 1],[0 1],[0 1]);
      grid.processGeometry;
      testCase.verifyEqual(grid.cells.num,5^3);
      testCase.verifyEqual(grid.nNodes,6^3);
      testCase.verifyEqual(grid.getCellNodes(2),[2	3	9	8	38	39	45	44]);

      b = BlockStructuredMesh(2,2,2,[0 1],[0 1],[0 1],3);
      b.refineRecursive([1 1 1],2);
      b.refineRecursive([2 2 2],2);
      b.refineRecursive([2 1 2],1);
      grid = b.processGeometry();
      testCase.verifyEqual(grid.nNodes,305);
      testCase.verifyEqual(grid.nDim,3);
      testCase.verifyEqual(grid.cells.num,141);
      testCase.verifyEqual(grid.getCellNodes(10),[28	32	33	29	30	34	35	31]);
    end
   end
  

end
