classdef IntegratedTests < matlab.unittest.TestCase

  properties
    testPath
    tWall
    tCPU
    oldVerbosity
  end

  methods (TestClassSetup)
    % Shared setup for the entire test class
    function startupGReS(testCase)
      % Make initGReS visible to this worker
      gresRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
      addpath(gresRoot);
      initGReS(0);

      testCase.testPath = fullfile(gres_root(), 'Code', 'integratedTest');
      addpath(testCase.testPath);

      % Store current verbosity and disable output during tests
      testCase.oldVerbosity = gresLog().getVerbosity();
      gresLog().setVerbosity(-2);
    end
  end

  methods (TestClassTeardown)
    % Executed once after all tests
    function shutdownGReS(testCase)
      gresLog().setVerbosity(testCase.oldVerbosity);
    end
  end

  methods (TestMethodSetup)
    % Setup executed before every test
    function setupTest(testCase)
      % cd(testCase.testPath);
      testCase.tWall = tic();
      testCase.tCPU  = cputime;
    end
  end

  methods (TestMethodTeardown)
    % Teardown executed after every test
    function teardownTest(testCase)
      elapsedWall = toc(testCase.tWall);
      elapsedCPU  = cputime - testCase.tCPU;

      fprintf( 'Wall time: %.3f s | CPU time: %.3f s\n', ...
        elapsedWall, ...
        elapsedCPU);

      close all;
      % cd(testCase.testPath);
    end
  end

  methods (Test)
    function Terzaghi(testCase)
      testFile = fullfile(testCase.testPath,'Terzaghi','testTerzaghi.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'Terzaghi tests did not pass.');
    end

    function SubDomains(testCase)
      testFile = fullfile(testCase.testPath,'SubDomains','testSubDomains.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'SubDomains tests did not pass.');
    end

    function Richards(testCase)
      testFile = fullfile(testCase.testPath,'Richards','testRichards.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'Richards tests did not pass.');
    end

    function CubeSedimentation(testCase)
      testFile = fullfile(testCase.testPath,'CubeSedimentation','testSedimentation.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'CubeSedimentation tests did not pass.');
    end

    function MortarConvergence(testCase)
      testFile = fullfile(testCase.testPath,'MortarConvergence','testMortarPoisson.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'MortarConvergence tests did not pass.');
    end

    function ConstantSliding(testCase)
      testFile = fullfile(testCase.testPath,'ConstantSliding','testConstantSliding.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'ConstantSliding tests did not pass.');
    end

    function SingleCrackCompressed(testCase)
      testFile = fullfile(testCase.testPath,'SingleCrackCompressed','testSingleCrackCompressed.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'SingleCrackCompressed tests did not pass.');
    end

    function ConstantSlidingEFEM(testCase)
      testFile = fullfile(testCase.testPath,'ConstantSlidingEFEM','testConstantSlidingEFEM.m');
      results = runtests(testFile);
      testCase.verifyFalse(any([results.Failed]),'ConstantSlidingEFEM tests did not pass.');
    end

  end
end