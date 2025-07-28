classdef testMaterialN < matlab.unittest.TestCase

   methods (TestClassSetup)
      % Shared setup for the entire test class
   end

   methods (TestMethodSetup)
      % Setup for each test
   end

   methods (Test)
      % Test methods

      function MaterialNTest(testCase)

         model = [];
         fdata = 'Inputs/materials.xml';
         m = MaterialN(model,fdata);             % Create material database

         pt=1;

         % testCase.verifyEqual(m.getData('steel'),struct('E', 200e9, 'rho', 7850))
      end
   end

end