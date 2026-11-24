% Run the input-free MortarConvergence study after initializing GReS.
assertSuccess(runtests(fullfile(fileparts(mfilename('fullpath')),'testMortarPoisson.m')));
