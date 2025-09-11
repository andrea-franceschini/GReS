classdef CheckTests < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods
        function Deep_Aquifer(testCase)
            cd("Deep_Aquifer/");
            run('Main.m');
            Sol_Press = [domain.outstate.results.expPress];
            Sol_Displ = [domain.outstate.results.expDispl];
            load("expData.mat")
            testCase.verifyEqual(Sol_Press,expPress);
            testCase.verifyEqual(Sol_Displ,expDispl);
            close all;
            cd("../");
        end

        function Mandel_Biot(testCase)
            cd("Mandel_Biot/");
            run('Main.m');
            Sol_Press = [domain.outstate.results.expPress];
            Sol_Displ = [domain.outstate.results.expDispl];
            load("expData.mat")
            testCase.verifyEqual(Sol_Press,expPress);
            testCase.verifyEqual(Sol_Displ,expDispl);
            close all;
            cd("../");
        end

        function Richards_Case1(testCase)
            cd("Richards/Case1/");
            run('Main.m');
            load("Inputs/Solution/UnitTest.mat")
            testCase.verifyEqual(pressplot,PressureUnity);
            testCase.verifyEqual(swplot,SaturationUnity);
            close all;
            cd("../../");
        end

        function Richards_Case2(testCase)
            cd("Richards/Case2/");
            run('Main.m');
            load("Inputs/Solution/UnitTest.mat")
            testCase.verifyEqual(pressplot,PressureUnity);
            testCase.verifyEqual(swplot,SaturationUnity);
            close all;
            cd("../../");
        end

        function Terzaghi_Biot(testCase)
            cd("Terzaghi_Biot/");
            run('Main.m');
            Sol_Press = printUtils.results.expPress;
            Sol_Displ = printUtils.results.expDispl;
            load("expData.mat")
            testCase.verifyEqual(Sol_Press,expPress);
            testCase.verifyEqual(Sol_Displ,expDispl);            
            close all;
            cd("../");
        end
    end

end