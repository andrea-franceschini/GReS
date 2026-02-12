classdef testMesh < matlab.unittest.TestCase

   methods(Test)
    % Test methods
    function testReadVTK(testCase)
      mesh = Mesh();
      mesh.importMesh('testMesh.vtk');
      testCase.verifyEqual(mesh.nNodes,919);
      testCase.verifyEqual(mesh.nDim,3);
      testCase.verifyEqual(mesh.nCells,560);
      testCase.verifyEqual(mesh.nSurfaces,624);
    end

    function testStructuredGrids(testCase)
      mesh = structuredMesh(5,5,5,[0 1],[0 1],[0 1]);
      testCase.verifyEqual(mesh.nCells,5^3);
      testCase.verifyEqual(mesh.nNodes,6^3);
      testCase.verifyEqual(mesh.cells(2,:),[2	3	9	8	38	39	45	44]);

      b = BlockStructuredMesh(2,2,2,[0 1],[0 1],[0 1],3);
      b.refineRecursive([1 1 1],2);
      b.refineRecursive([2 2 2],2);
      b.refineRecursive([2 1 2],1);
      mesh =b.processGeometry();
      testCase.verifyEqual(mesh.nNodes,305);
      testCase.verifyEqual(mesh.nDim,3);
      testCase.verifyEqual(mesh.nCells,141);
      testCase.verifyEqual(mesh.cells(10,:),[28	32	33	29	30	34	35	31]);
    end
   end
  

end