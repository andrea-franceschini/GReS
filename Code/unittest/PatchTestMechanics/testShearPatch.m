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
      gresLog().setVerbosity(-2);
      simparams = SimulationParameters(testCase.pathToFile);
      mesh = structuredMesh(1,1,1,[0 1],[0 1],[0 1]);
      elems = Elements(mesh,2);
      faces = Faces(mesh);
      grid = struct('topology',mesh,'cells',elems,'faces',faces);
      mat = Materials(testCase.pathToFile);
      bc = Boundaries(testCase.pathToFile,grid);
      domain = Discretizer('Boundaries',bc,...
                            'Materials',mat,...
                            'Grid',grid);
      domain.addPhysicsSolver(testCase.pathToFile);
      solver = NonLinearImplicit('simulationparameters',simparams,'domain',domain);
      solver.simulationLoop();
      gresLog().setVerbosity(-1);
      verifyEqual(testCase,domain.state.data.stress(:,5),1.5*ones(8,1),"AbsTol",1e-9)


    end
  end

end