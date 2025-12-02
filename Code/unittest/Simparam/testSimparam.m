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
  end

end