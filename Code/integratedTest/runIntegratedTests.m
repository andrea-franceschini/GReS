function results = runIntegratedTests(useParallel)
    % Run GReS integrated tests with optional parallel execution.
    % Usage:
    %   runIntegratedTests()       % Runs in parallel by default
    %   runIntegratedTests(false)  % Runs sequentially
    %   runIntegratedTests(true)   % Runs in parallel

    if nargin < 1
        useParallel = true;
    end

<<<<<<< HEAD
    testPath = fileparts(mfilename('fullpath'));
    testFile = fullfile(testPath, 'IntegratedTests.m');
=======
v = gresLog().getVerbosity();
gresLog().setVerbosity(-2);
>>>>>>> a89efce9a482544a17704feb06040edeee9becaa

    % Optional: manage parallel pool lifecycle if needed
    % if useParallel && isempty(gcp('nocreate'))
    %     parpool('Processes', 2);
    % end

    results = runtests(testFile, 'UseParallel', useParallel);

    disp(results);

<<<<<<< HEAD
    if any([results.Failed])
        error('GReS:IntegratedTestsFailed', 'One or more integrated tests failed.');
    end
end
=======
fprintf("Elapsed wall-clock time: %1.2f s\n", elapsedWall);
fprintf("Elapsed CPU time:        %1.2f s\n", elapsedCPU);

if any([results.Failed])
    error("Some tests did not pass");
else
    disp("All tests ran successfully");
end

gresLog().setVerbosity(v);
>>>>>>> a89efce9a482544a17704feb06040edeee9becaa
