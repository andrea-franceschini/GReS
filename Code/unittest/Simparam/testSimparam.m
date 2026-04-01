classdef testSimparam < matlab.unittest.TestCase

  properties
    model
    simParam
  end

  methods (TestClassSetup)
    % Shared setup for the entire test class
  end

  methods (TestMethodSetup)
    % Setup for each test
  end

  methods (Test)
    % Test methods

    function SimparamXmlTest(testCase)
      fileName = 'Input/simParam.xml';
      testCase.simParam = SimulationParameters(fileName);
      verifyEqual(testCase,testCase.simParam.tMax,1e1);
      verifyEqual(testCase,testCase.simParam.itMaxNR,10);
      verifyEqual(testCase,testCase.simParam.relTol,1e-8);
      verifyEqual(testCase,testCase.simParam.absTol,1e-9);
    end

    function SimparamKVTest(testCase)
      testCase.simParam = SimulationParameters('Start',0.0,...
                                               'End',1e1,...
                                               'DtInit',0.5,...
                                               'DtMax',5.0,...
                                               'DtMin',0.5,...
                                               'incrementFactor',1.1,...
                                               'MaxNLIteration',10,....
                                               'RelativeTolerance',1.e-8,...
                                               'AbsoluteTolerance',1.e-9);
      verifyEqual(testCase,testCase.simParam.tMax,1e1);
      verifyEqual(testCase,testCase.simParam.itMaxNR,10);
      verifyEqual(testCase,testCase.simParam.relTol,1e-8);
      verifyEqual(testCase,testCase.simParam.absTol,1e-9);
    end
  end

end