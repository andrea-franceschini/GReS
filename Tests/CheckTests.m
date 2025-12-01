classdef CheckTests < matlab.unittest.TestCase

  methods (TestClassSetup)
    %Shared setup for the entire test class
    function startupGReS(testCase)
      addpath('../')
      initGReS(0);                   % main MATLAB
    end
  end

  methods (TestClassSetup)

  end

  methods (Test)
    % Test methods
    function Deep_Aquifer(testCase)
      cd("Deep_Aquifer/");
      run('Main.m');
      % checking only if 3 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),3);
      testCase.verifyEqual(domain.outstate.timeList,[1; 4; 8]);
      close all;
      cd("../");
    end

    function Mandel_Biot(testCase)
      cd("Mandel_Biot/");
      run('Main.m');
      % checking only if 5 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),5);
      testCase.verifyEqual(domain.outstate.timeList,[0.0500; 0.2500; 1; 2.5000; 5]);
      close all;
      cd("../");
    end

    function Richards_Case1(testCase)
      cd("Richards/Case1/");
      run('Main.m');
      % checking only if 3 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),3);
      testCase.verifyEqual(domain.outstate.timeList,[10; 50; 100]);
      close all;
      cd("../../");
    end

    function Richards_Case2(testCase)
      cd("Richards/Case2/");
      run('Main.m');
      % checking only if 4 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),4);
      testCase.verifyEqual(domain.outstate.timeList,[51840; 77760; 129600; 259200]);
      close all;
      cd("../../");
    end

    function Terzaghi_Biot(testCase)
      cd("Terzaghi_Biot/");
      run('Main.m');
      % checking only if 6 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),6);
      testCase.verifyEqual(domain.outstate.timeList,[15; 30; 60; 90; 120; 180]);
      close all;
      cd("../");
    end

    function Flux_barrier(testCase)
      cd("Flux_barrier/");
      run('Main.m');
      % checking only if 2 outputs is being generated
      testCase.verifyEqual(length(domain.outstate.results),1);
      testCase.verifyEqual(domain.outstate.timeList,10);
      close all;
      cd("../");
    end

    function Flow_nonConforming(testCase)
      cd("Flow_nonConforming/");
      main
      % checking only if 6 outputs is being generated
      testCase.verifyEqual(domains(1).outstate.timeList,domains(2).outstate.timeList);
      testCase.verifyEqual(domains(1).outstate.timeList,[0.5000; 1; 2; 5; 10; 50]);
      close all;
      cd("../");
    end
  end

end