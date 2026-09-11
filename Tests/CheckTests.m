classdef CheckTests < matlab.unittest.TestCase

  properties
    testPath
    tWall
    tCPU
  end

  methods (TestClassSetup)
    % Shared setup for the entire test class
    function startupGReS(testCase)
      % Make initGReS visible to this worker
      addpath(fileparts(fileparts(mfilename('fullpath'))));
      initGReS(0);
      testCase.testPath = fullfile(gres_root(), 'Tests');
      addpath(testCase.testPath);
      gresLog().setVerbosity(-2);
    end
  end

  methods (TestMethodSetup)
    % Setup executed before every test
    function setupTest(testCase)
      cd(testCase.testPath);
      testCase.tWall = tic();
      testCase.tCPU  = cputime;
    end
  end

  methods (TestMethodTeardown)
    % Teardown executed after every test.
    function teardownTest(testCase)
      elapsedWall = toc(testCase.tWall);
      elapsedCPU  = cputime - testCase.tCPU;
      fprintf('Wall time: %.3f s | CPU time: %.3f s\n',elapsedWall,elapsedCPU);
      close all;
      cd(testCase.testPath);
    end
  end

  methods (Test)
    % Test methods
    function DeepAquifer(testCase)
      cd("DeepAquifer/");
      run('runDeepAquifer.m');
      % checking only if 3 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),6);
      % testCase.verifyEqual(domain.outstate.timeList,[1; 4; 8]);
    end

    function MandelBiot(testCase)
      cd("MandelBiot/");
      run('runMandelBiot.m');
      % checking only if 5 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),5);
      % testCase.verifyEqual(domain.outstate.timeList',[0.0500; 0.2500; 1; 2.5000; 5]);
    end

    function RichardsCase1(testCase)
      cd("Richards/Case1/");
      run('runRichardCase1.m');
      % checking only if 3 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),3);
      testCase.verifyEqual(domain.outstate.timeList',[10; 50; 100]);
    end

    function RichardsCase2(testCase)
      cd("Richards/Case2/");
      run('runRichardCase2.m');
      % checking only if 4 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),4);
      testCase.verifyEqual(domain.outstate.timeList',[51840; 77760; 129600; 259200]);
    end

    function TerzaghiBiot(testCase)
      cd("TerzaghiBiot/");
      run('runTerzaghi.m');
      % checking only if 6 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),6);
      % testCase.verifyEqual(domain.outstate.timeList,[15; 30; 60; 90; 120; 180]);
    end

    function FluxBarrier(testCase)
      cd("FluxBarrier/");
      run('runFluxBarrier.m');
      % checking only if 2 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),1);
      % testCase.verifyEqual(domain.outstate.timeList,10);
    end

    function FlowNonConforming(testCase)
      cd("FlowNonConforming/");
      run('runFlowNonConforming.m');
      % checking only if 6 outputs is being generated
      % testCase.verifyEqual(domains(1).outstate.timeList,domains(2).outstate.timeList);
      % testCase.verifyEqual(domains(1).outstate.timeList,[0.5000; 1; 2; 5; 10; 50]);
    end

    function StickSlipOpen(testCase)
      cd("StickSlipOpenContact/");
      run('runStickSlipOpen.m');
    end

    function SneddonEFEM(testCase)
      cd("SneddonProblemEFEM/");
      run('runSneddon.m');
    end
  end

end
