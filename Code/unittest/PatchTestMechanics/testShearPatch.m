classdef testShearPatch < matlab.unittest.TestCase
  
  properties
    pathToFile
  end

  methods(TestClassSetup)
    % Shared setup for the entire test class
    function setupOnce(testCase)
      testCase.pathToFile = fullfile(gres_root,...
        "/Code/unittest/PatchTestMechanics/shearPatch.xml");
    end
  end
  
  methods(TestMethodSetup)
    % Setup for each test
  end

  methods(Test)
    % Test methods

    function testPatch(testCase)
      gresLog().setVerbosity(-1);
      simparams = SimulationParameters(testCase.pathToFile);
      b = BlockStructuredMesh([0.0, 1.0;0.0 1.0; 0.0 1.0],[1,1,1],1);
      mesh = b.processGeometry();
      elems = Elements(mesh,2);
      faces = Faces(mesh);
      grid = struct('topology',mesh,'cells',elems,'faces',faces);
      mat = Materials(testCase.pathToFile);
      printUtils = OutState(mesh,"writeVtk",0,"flagMatFile",0);
      bc = Boundaries(testCase.pathToFile,grid);
      domain = Discretizer('Boundaries',bc,...
        'OutState',printUtils,...
        'Materials',mat,...
        'Grid',grid);
      domain.addPhysicsSolver(testCase.pathToFile);
      solver = GeneralSolver(simparams,domain);
      solver.NonLinearLoop();

      verifyEqual(testCase,domain.state.data.stress(:,5),1.5*ones(8,1),"AbsTol",1e-9)


    end
  end

end