classdef testBoundaries < matlab.unittest.TestCase
  
  properties
    bc
  end

  methods(TestClassSetup)
    % Shared setup for the entire test class
    function setupOnce(testCase)
      testCase.pathToFile = fullfile(gres_root,"/Code/unittest/BoundaryConditions");
    end
  end
  
  methods(TestMethodSetup)
    % Setup for each test
  end
  
  methods(Test)
    % Test methods
    
    function readBC(testCase)
      model = ModelType(["Poromechanics_FEM","SinglePhaseFlow_FVTPFA"]);
      mesh.importMesh("Column_hexa.msh");
      elems = Elements(mesh,2);
      faces = Faces(model,mesh);
      grid = struct("topology",mesh,"cells",elems,"faces",faces);
      testCase.bc = BoundariesNew(testCase.pathToFile,model,grid);
      
      v1 = testCase.bc.getValue("bc1");
      v2 = testCase.bc.getValue("bc2");
      v3 = testCase.bc.getValue("bc3");
      v4 = testCase.bc.getValue("bc4");
      v5 = testCase.bc.getValue("bc5");

      e1 = testCase.bc.getEntities("bc1");
      e2 = testCase.bc.getEntities("bc2");
      e3 = testCase.bc.getEntities("bc3");
      e4 = testCase.bc.getEntities("bc4");
      e5 = testCase.bc.getEntities("bc5");

      le2 = testCase.bc.getLoadedEnts("bc2"); 



     
    end
  end
  
end