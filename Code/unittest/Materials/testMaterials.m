classdef testMaterials < matlab.unittest.TestCase

  properties
    mat
    pathToFile
  end

  methods (TestClassSetup)
    % Shared setup for the entire test class
    function setupOnce(testCase)
      testCase.pathToFile = fullfile(gres_root,"Code/unittest/Materials/Input/materials.xml");
    end
  end

  methods (TestMethodSetup)
    % Setup for each test
  end

  methods (Test)
    % Test methods

    function readFile(testCase)
      testCase.mat = MaterialsXML(testCase.pathToFile);
      verifyEqual(testCase,getMaterial(testCase.mat,1).constLaw.E,1e3);
      K = diag([1e-12;1e-12;1e-12]);
      verifyEqual(testCase,getMaterial(testCase.mat,4).PorousRock.getPermMatrix(),K);
      verifyEqual(testCase,getMaterial(testCase.mat,2).constLaw.phi,30);
      verifyEqual(testCase,getFluid(testCase.mat).getFluidSpecWeight(),0);
      verifyEqual(testCase,getFluid(testCase.mat).getFluidCompressibility(),4.4e-7);
    end
  end

end