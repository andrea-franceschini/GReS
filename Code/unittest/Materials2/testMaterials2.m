classdef testMaterials2 < matlab.unittest.TestCase

  properties
    mat
    pathToFile
  end

  methods (TestClassSetup)
    % Shared setup for the entire test class
    function setupOnce(testCase)
      testCase.pathToFile = fullfile(gres_root,"Code/unittest/Materials2/Input/materials.xml");
    end
  end

  methods (TestMethodSetup)
    % Setup for each test
  end

  methods (Test)
    % Test methods

    function readFile(testCase)
        grid = structuredMesh(5,5,5,[0 1],[0 1],[0 1]);
        grid.processGeometry;
        mat = MaterialS(grid);
        mat.registerParameter('permeability',ProprietyType.scalar,);
    end

    % function readFile(testCase)
    %   testCase.mat = Materials(testCase.pathToFile);
    %   verifyEqual(testCase,getMaterial(testCase.mat,1).ConstLaw.E,1e3);
    %   K = diag([1e-12;1e-12;1e-12]);
    %   verifyEqual(testCase,getMaterial(testCase.mat,4).PorousRock.getPermMatrix(),K);
    %   verifyEqual(testCase,getMaterial(testCase.mat,2).ConstLaw.phi,30);
    %   verifyEqual(testCase,getFluid(testCase.mat).getSpecificWeight(),0);
    %   verifyEqual(testCase,getFluid(testCase.mat).getFluidCompressibility(),4.4e-7);
    % end
    % 
    % function setMaterial(testCase)
    %   testCase.mat = Materials();
    %   testCase.mat.addSolid('name',"rock",'cellTags',1);
    %   testCase.mat.addConstitutiveLaw("rock","Elastic",'youngModulus',5e3,'poissonRatio',0.25);
    %   testCase.mat.addFluid('dynamicViscosity',1e-3);
    %   testCase.mat.addPorousRock("rock","specificWeight",21.0,"permeability",1e-12);
    %   testCase.mat.addCapillaryCurves("rock","type","mualem","beta",2.0,"n",1.0,"kappa",1.0)
    %   verifyEqual(testCase,getMaterial(testCase.mat,1).ConstLaw.E,5e3);
    %   K = diag([1e-12;1e-12;1e-12]);
    %   verifyEqual(testCase,getMaterial(testCase.mat,1).PorousRock.getPermMatrix(),K);
    %   verifyEqual(testCase,getMaterial(testCase.mat,1).PorousRock.Curves.beta,2.0);
    % end
    % 
    % function testStructuredGrids(testCase)
    %   grid = structuredMesh(5,5,5,[0 1],[0 1],[0 1]);
    %   grid.processGeometry;
    %   testCase.verifyEqual(grid.cells.num,5^3);
    %   testCase.verifyEqual(grid.nNodes,6^3);
    %   testCase.verifyEqual(grid.getCellNodes(2),[2	3	9	8	38	39	45	44]);
    % 
    %   b = BlockStructuredMesh(2,2,2,[0 1],[0 1],[0 1],3);
    %   b.refineRecursive([1 1 1],2);
    %   b.refineRecursive([2 2 2],2);
    %   b.refineRecursive([2 1 2],1);
    %   grid = b.processGeometry();
    %   testCase.verifyEqual(grid.nNodes,305);
    %   testCase.verifyEqual(grid.nDim,3);
    %   testCase.verifyEqual(grid.cells.num,141);
    %   testCase.verifyEqual(grid.getCellNodes(10),[28	32	33	29	30	34	35	31]);
    % end



  end

end