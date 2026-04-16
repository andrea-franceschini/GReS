classdef testMesh < matlab.unittest.TestCase

   methods(Test)
    % Test methods
    function testReadVTK(testCase)
      grid = Grid();
      grid.importMesh('testMesh.vtk');
      testCase.verifyEqual(grid.nNodes,919);
      testCase.verifyEqual(grid.nDim,3);
      testCase.verifyEqual(grid.cells.num,560);
      testCase.verifyEqual(grid.surfaces.num,624);
    end

    function testStructuredGrids(testCase)
      grid = structuredMesh(5,5,5,[0 1],[0 1],[0 1]);
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