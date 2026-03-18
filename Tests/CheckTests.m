classdef CheckTests < matlab.unittest.TestCase

  methods (TestClassSetup)
    %Shared setup for the entire test class
    function startupGReS(testCase)
      addpath('../')
      initGReS(0);                   % main MATLAB
      gresLog().setVerbosity(-2)
    end
  end

  methods (TestClassSetup)

  end

  methods (Test)
    % Test methods
    function DeepAquifer(testCase)
      cd("DeepAquifer/");
      run('Main.m');
      % checking only if 3 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),6);
      % testCase.verifyEqual(domain.outstate.timeList,[1; 4; 8]);
      close all;
      cd("../");
    end

    function MandelBiot(testCase)
      cd("MandelBiot/");
      run('Main.m');
      % checking only if 5 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),5);
      % testCase.verifyEqual(domain.outstate.timeList',[0.0500; 0.2500; 1; 2.5000; 5]);
      close all;
      cd("../");
    end

    function RichardsCase1(testCase)
      cd("Richards/Case1/");
      run('Main.m');
      % checking only if 3 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),3);
      testCase.verifyEqual(domain.outstate.timeList',[10; 50; 100]);
      close all;
      cd("../../");
    end

    function RichardsCase2(testCase)
      cd("Richards/Case2/");
      run('Main.m');
      % checking only if 4 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),4);
      testCase.verifyEqual(domain.outstate.timeList',[51840; 77760; 129600; 259200]);
      close all;
      cd("../../");
    end

    function TerzaghiBiot(testCase)
      cd("TerzaghiBiot/");
      run('Main.m');
      % checking only if 6 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),6);
      % testCase.verifyEqual(domain.outstate.timeList,[15; 30; 60; 90; 120; 180]);
      close all;
      cd("../");
    end

    function FluxBarrier(testCase)
      cd("FluxBarrier/");
      run('Main.m');
      % checking only if 2 outputs is being generated
      % testCase.verifyEqual(length(domain.outstate.results),1);
      % testCase.verifyEqual(domain.outstate.timeList,10);
      close all;
      cd("../");
    end

    function FlowNonConforming(testCase)
      cd("FlowNonConforming/");
      run('Main.m');
      % checking only if 6 outputs is being generated
      % testCase.verifyEqual(domains(1).outstate.timeList,domains(2).outstate.timeList);
      % testCase.verifyEqual(domains(1).outstate.timeList,[0.5000; 1; 2; 5; 10; 50]);
      close all;
      cd("../");
    end

    function StickSlipOpen(testCase)
      cd("StickSlipOpenContact/");
      run('StickSlipOpen.m');
      close all;
      cd("../");
    end

    function SneddonEFEM(testCase)
      cd("SneddonProblemEFEM/");
      run('sneddon.m');
      close all;
      cd("../");
    end
  end

end