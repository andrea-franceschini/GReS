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
   end
  

end