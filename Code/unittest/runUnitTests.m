% Run GReS unit tests in parallel.
testPath = fullfile(fileparts(mfilename('fullpath')));

% pool = gcp('nocreate');
% if isempty(pool)
%     parpool('Processes', 2);
% end

results = runtests(fullfile(testPath, 'UnitTests.m'),'UseParallel', true);
disp(results);
if any([results.Failed])
  error('GReS:UnitTestsFailed','One or more unit tests failed.');
else
  disp('All unit tests passed.');
end