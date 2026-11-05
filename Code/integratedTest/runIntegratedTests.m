function results = runIntegratedTests(useParallel)
    % Run GReS integrated tests with optional parallel execution.
    % Usage:
    %   runIntegratedTests()       % Runs in parallel by default
    %   runIntegratedTests(false)  % Runs sequentially
    %   runIntegratedTests(true)   % Runs in parallel

    if nargin < 1
        useParallel = true;
    end

    testPath = fileparts(mfilename('fullpath'));
    testFile = fullfile(testPath, 'IntegratedTests.m');

    % Optional: manage parallel pool lifecycle if needed
    % if useParallel && isempty(gcp('nocreate'))
    %     parpool('Processes', 2);
    % end

    results = runtests(testFile, 'UseParallel', useParallel);

    disp(results);

    if any([results.Failed])
        error('GReS:IntegratedTestsFailed', 'One or more integrated tests failed.');
    end
end