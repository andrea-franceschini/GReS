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
      testCase.mat = Materials(testCase.pathToFile);
      verifyEqual(testCase,getMaterial(testCase.mat,1).ConstLaw.E,1e3);
      K = diag([1e-12;1e-12;1e-12]);
      verifyEqual(testCase,getMaterial(testCase.mat,4).PorousRock.getPermMatrix(),K);
      verifyEqual(testCase,getMaterial(testCase.mat,2).ConstLaw.phi,30);
      verifyEqual(testCase,getFluid(testCase.mat).getSpecificWeight(),0);
      verifyEqual(testCase,getFluid(testCase.mat).getFluidCompressibility(),4.4e-7);
    end

    function setMaterial(testCase)
      testCase.mat = Materials();
      testCase.mat.addSolid('name',"rock",'cellTags',1);
      testCase.mat.addConstitutiveLaw("rock","Elastic",'youngModulus',5e3,'poissonRatio',0.25);
      testCase.mat.addFluid('dynamicViscosity',1e-3);
      testCase.mat.addPorousRock("rock","specificWeight",21.0,"permeability",1e-12);
      testCase.mat.addCapillaryCurves("rock","type","mualem","beta",2.0,"n",1.0,"kappa",1.0)
      verifyEqual(testCase,getMaterial(testCase.mat,1).ConstLaw.E,5e3);
      K = diag([1e-12;1e-12;1e-12]);
      verifyEqual(testCase,getMaterial(testCase.mat,1).PorousRock.getPermMatrix(),K);
      verifyEqual(testCase,getMaterial(testCase.mat,1).PorousRock.Curves.beta,2.0);
    end
  end

end