% Run GReS integrated tests in parallel.
testPath = fullfile(fileparts(mfilename('fullpath')));

% pool = gcp('nocreate');
% if isempty(pool)
%     parpool('Processes', 2);
% end

results = runtests(fullfile(testPath, 'IntegratedTests.m'),'UseParallel', true);

disp(results);

if any([results.Failed])
    error('GReS:IntegratedTestsFailed','One or more integrated tests failed.');
end