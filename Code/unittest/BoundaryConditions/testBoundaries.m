classdef testBoundaries < matlab.unittest.TestCase
  
  properties
    bc
    pathToFile
  end

  methods(TestClassSetup)
    % Shared setup for the entire test class
    function setupOnce(testCase)
      testCase.pathToFile = fullfile(gres_root,...
        "/Code/unittest/BoundaryConditions/boundaryConditions.xml");
    end
  end
  
  methods(TestMethodSetup)
    % Setup for each test
  end
  
  methods(Test)
    % Test methods
    
    function readBC(testCase)
      model = ModelType(["Poromechanics_FEM","SinglePhaseFlow_FVTPFA"]);
      mesh = Mesh();
      mesh.importMesh("mesh.msh");
      elems = Elements(mesh,2);
      faces = Faces(model,mesh);
      grid = struct("topology",mesh,"cells",elems,"faces",faces);
      testCase.bc = Boundaries(testCase.pathToFile,model,grid);

      v1 = testCase.bc.getVals("bc1",3);
      v21 = testCase.bc.getVals("bc2",0.2);
      v22 = testCase.bc.getVals("bc2",2);
      v23 = testCase.bc.getVals("bc2",10);
      v3 = testCase.bc.getVals("bc3",1.1);
      v4 = testCase.bc.getVals("bc4",2);
      v5 = testCase.bc.getVals("bc5",0.0);


      verifyEqual(testCase,v1,[60;34.50;0.0],"AbsTol",1e-8)
      verifyEqual(testCase,mean(v21),-10.0,"AbsTol",1e-9)
      verifyEqual(testCase,mean(v22),-7.5,"AbsTol",1e-9)
      verifyEqual(testCase,mean(v23),10.0,"AbsTol",1e-9)
      verifyEqual(testCase,length(v3),126,"AbsTol",1e-9)
      verifyEqual(testCase,mean(v3),5,"AbsTol",1e-9)
      verifyEqual(testCase,length(v4),126,"AbsTol",1e-9)
      verifyEqual(testCase,max(v5),0,"AbsTol",1e-9)
      verifyEqual(testCase,length(v5),12,"AbsTol",1e-9)


      e1 = testCase.bc.getEntities("bc1");
      e2 = testCase.bc.getEntities("bc2");
      e3 = testCase.bc.getEntities("bc3");
      e5 = testCase.bc.getEntities("bc5");

      verifyEqual(testCase,e1,[1;54;20],"AbsTol",1e-9)
      verifyEqual(testCase,e2,repmat((165:168)',2,1),"AbsTol",1e-9)
      verifyEqual(testCase,length(e3),126,"AbsTol",1e-9)
      verifyEqual(testCase,e5,repmat((1:4)',3,1),"AbsTol",1e-9)


      le2 = testCase.bc.getLoadedEntities("bc2"); 
      leCheck = repmat([5;6;7;8;13;14;15;16;170],2,1);
      verifyEqual(testCase,le2,leCheck,"AbsTol",1e-9);

      loadArea = testCase.bc.getEntitiesInfluence("bc2");
      area = nonzeros(loadArea);
      verifyEqual(testCase,area,repmat(0.0625,32,1),"AbsTol",1e-9);

    end
  end
  
end