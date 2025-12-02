suite = matlab.unittest.TestSuite.fromClass(?CheckTests);
runner = matlab.unittest.TestRunner.withTextOutput;
results = runner.runInParallel(suite);
disp(results)