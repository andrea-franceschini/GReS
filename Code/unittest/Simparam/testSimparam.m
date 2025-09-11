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
    function SimparamDatTest1(testCase)
      fileName = 'Input/SimParamSPF.dat';
      testCase.model = ModelType("SinglePhaseFlow_FVTPFA");
      testCase.simParam = SimulationParameters(fileName,testCase.model);
      verifyEqual(testCase,testCase.simParam.tMax,100);
      verifyEqual(testCase,testCase.simParam.itMaxNR,12);
      verifyEqual(testCase,testCase.simParam.relTol,1e-6);
      verifyEqual(testCase,testCase.simParam.absTol,1e-10);
    end

    function SimparamDatTest2(testCase)
      fileName = 'Input/SimParamVSF.dat';
      testCase.model = ModelType("VariabSatFlow_FVTPFA");
      testCase.simParam = SimulationParameters(fileName,testCase.model);
      verifyEqual(testCase,testCase.simParam.tMax,100);
      verifyEqual(testCase,testCase.simParam.itMaxNR,10);
      verifyEqual(testCase,testCase.simParam.NLSolver,'Picard');
      verifyEqual(testCase,testCase.simParam.relTol,1e-6);
      verifyEqual(testCase,testCase.simParam.absTol,1e-10);
    end

    function SimparamXmlTest(testCase)
      fileName = 'Input/SimParam.xml';
      testCase.model = ModelType("VariabSatFlow_FVTPFA");
      testCase.simParam = SimulationParameters(fileName,testCase.model);
      verifyEqual(testCase,testCase.simParam.tMax,259200);
      verifyEqual(testCase,testCase.simParam.itMaxNR,10);
      verifyEqual(testCase,testCase.simParam.NLSolver,'Newton');
      verifyEqual(testCase,testCase.simParam.relTol,1e-8);
      verifyEqual(testCase,testCase.simParam.absTol,1e-9);
    end
  end

end