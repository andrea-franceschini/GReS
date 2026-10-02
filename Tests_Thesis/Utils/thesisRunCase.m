function report = thesisRunCase(test, options)
% Never add all case directories with genpath: helper names overlap.
caseRoot = fullfile(options.SuiteRoot, sprintf('Chapter_%02d', test.chapter), test.name);
runDir = fullfile(caseRoot, 'Output');
if isfolder(runDir)
    rmdir(runDir, 's');
end
mkdir(runDir);
report = struct('name', string(test.name), 'section', string(test.section), ...
    'status', "pending", 'seconds', 0, 'output', string(runDir), 'message', "");
if ~isempty(test.blocked)
    report.status = "blocked";
    report.message = string(test.blocked);
    save(fullfile(runDir, 'status.mat'), 'report');
    return
end
oldDir = pwd;
oldPath = path;
oldStream = RandStream.getGlobalStream;
oldStreamState = oldStream.State;
oldVisibility = get(groot, 'defaultFigureVisible');
oldFigures = findall(groot, 'Type', 'figure');
cleanup = onCleanup(@() restoreSession(oldDir, oldPath, oldStream, oldStreamState, oldVisibility, oldFigures)); %#ok<NASGU>
ticID = tic;
try
    % Case helpers take precedence over suite helpers; empty Utils is optional.
    addpath(fullfile(options.SuiteRoot, 'Utils'), '-end');
    caseUtils = fullfile(caseRoot, 'Utils');
    if isfolder(caseUtils)
        addpath(caseUtils, '-begin');
    end
    addpath(caseRoot, '-begin');
    cd(caseRoot);
    RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', options.Seed));
    if ~options.Visible
        set(groot, 'defaultFigureVisible', 'off');
    end
    diary(fullfile(runDir, 'run.log'));
    fprintf('<strong> ============================================================== </strong>\n\n');
    fprintf('<strong> Case: %s; section %s </strong> \n\n', test.title, test.section);
    fprintf('<strong> ============================================================== </strong>\n\n');
    %fprintf('GReS: %s\n', options.GReSRoot);
    assert(exist('mxGetDerBasisAndDet', 'file') == 3, ...
        'Thesis:MissingMEX', 'Run compileAll in the target GReS checkout before the suite.');
    result = main(options);
    save(fullfile(runDir, 'result.mat'), 'result', 'options', 'test', '-v7.3');
    report.status = "completed"; % execution completed; not a claim of thesis validation
    if isfield(result, 'failedConfigurations') && result.failedConfigurations > 0
        report.status = "completed_with_failures";
        report.message = sprintf('%d configurations failed; see result.mat.', result.failedConfigurations);
    end
catch ME
    report.status = "failed";
    report.message = string(getReport(ME, 'extended', 'hyperlinks', 'off'));
    fprintf(2, '%s\n', report.message);
end
report.seconds = toc(ticID);
diary off;
save(fullfile(runDir, 'status.mat'), 'report');
end

function restoreSession(oldDir, oldPath, oldStream, oldStreamState, oldVisibility, oldFigures)
diary off;
newFigures = setdiff(findall(groot, 'Type', 'figure'), oldFigures);
close(newFigures);
cd(oldDir);
path(oldPath);
oldStream.State = oldStreamState;
RandStream.setGlobalStream(oldStream);
set(groot, 'defaultFigureVisible', oldVisibility);
end
