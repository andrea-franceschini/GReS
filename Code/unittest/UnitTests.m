classdef UnitTests < matlab.unittest.TestCase

  properties
    testPath
    tWall
    tCPU
    oldVerbosity
  end

  methods (TestClassSetup)
    function startupGReS(testCase)
      % Make initGReS visible to this worker
      testFilePath = fileparts(mfilename('fullpath'));
      gresRoot = fileparts(fileparts(testFilePath));
      addpath(gresRoot);

      initGReS(0);
      testCase.testPath = fullfile(gres_root(),'Code','unittest');
      addpath(testCase.testPath);

      % Store current verbosity and disable output during tests
      testCase.oldVerbosity = gresLog().getVerbosity();
      gresLog().setVerbosity(-2);
    end
  end

  methods (TestClassTeardown)
    function shutdownGReS(testCase)
      gresLog().setVerbosity(testCase.oldVerbosity);
    end
  end

  methods (TestMethodSetup)
    function setupTest(testCase)
      testCase.tWall = tic();
      testCase.tCPU  = cputime;
    end
  end

  methods (TestMethodTeardown)
    function teardownTest(testCase)
      elapsedWall = toc(testCase.tWall);
      elapsedCPU  = cputime - testCase.tCPU;
      fprintf('Wall time: %.3f s | CPU time: %.3f s\n', ...
        elapsedWall,elapsedCPU);
      % close all;
    end
  end

  methods (Test)
    function Mesh(testCase)
      testFile = fullfile(testCase.testPath,'Mesh','testMesh.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'Mesh tests did not pass.');
    end

    function Simparam(testCase)
      testFile = fullfile(testCase.testPath,'Simparam','testSimparam.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'Simparam tests did not pass.');
    end

    function Materials(testCase)
      testFile = fullfile(testCase.testPath,'Materials','testMaterials.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'Materials tests did not pass.');
    end

    function BoundaryConditions(testCase)
      testFile = fullfile(testCase.testPath,'BoundaryConditions','testBoundaries.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'BoundaryConditions tests did not pass.');
    end

    function PatchTestMechanics(testCase)
      testFile = fullfile(testCase.testPath,'PatchTestMechanics','testShearPatch.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'PatchTestMechanics tests did not pass.');
    end

  end
end