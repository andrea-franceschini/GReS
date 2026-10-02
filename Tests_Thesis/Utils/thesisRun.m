function report = thesisRun(suiteRoot, varargin)
% Run selected thesis cases in sequence after initGReS.
p = inputParser;
addParameter(p, 'Chapter', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Cases', strings(0), @(x) isstring(x) || iscellstr(x) || ischar(x));
addParameter(p, 'Smoke', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Seed', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ListOnly', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Visible', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
options = p.Results;
registry = jsondecode(fileread(fullfile(suiteRoot, 'Input', 'testRegistry.json')));
if ~isempty(options.Chapter)
    registry = registry(ismember([registry.chapter], options.Chapter));
end
if ~isempty(options.Cases)
    requested = string(options.Cases);
    assert(all(ismember(requested, string({registry.name}))), ...
        'Thesis:UnknownCase', 'Unknown case name for this selection. Use ListOnly.');
    registry = registry(ismember(string({registry.name}), requested));
end
assert(~isempty(registry), 'Thesis:EmptySelection', 'No cases selected.');
% Smoke mode excludes cases without a smoke configuration before execution.
if options.Smoke
    supported = [registry.smokeSupported];
    for k = find(~supported(:)).'
        fprintf('Skipped %s: smoke mode is not supported.\n', registry(k).name);
    end
    registry = registry(supported);
    if isempty(registry)
        if options.ListOnly
            report = struct2table(registry);
        else
            emptyResult = struct('name', "", 'section', "", 'status', "", ...
                'seconds', 0, 'output', "", 'message', "");
            report = struct2table(repmat(emptyResult, 0, 1));
        end
        fprintf('No smoke-supported cases selected. Nothing was run.\n');
        return
    end
end
if options.ListOnly
    report = struct2table(registry);
    disp(report);
    return
end

% initGReS stores the selected checkout here; no path argument is needed.
gresRoot = getappdata(0, 'gres_root');
assert((ischar(gresRoot) || isstring(gresRoot)) && ...
    isfile(fullfile(gresRoot, 'initGReS.m')), ...
    'Thesis:GReSNotInitialized', 'Run initGReS once before running thesis tests.');
options.GReSRoot = gresRoot;
options.SuiteRoot = suiteRoot;

% The summary is replaced on every suite run.
summaryDir = fullfile(suiteRoot, 'Output');
if isfolder(summaryDir)
    rmdir(summaryDir, 's');
end
mkdir(summaryDir);
results = cell(numel(registry), 1);
for k = 1:numel(registry)
    results{k} = thesisRunCase(registry(k), options);
    fprintf('[%d/%d] %s: %s\n', k, numel(registry), results{k}.name, results{k}.status);
end
report = struct2table(vertcat(results{:}));
save(fullfile(summaryDir, 'report.mat'), 'report', 'options', 'registry');
writetable(report, fullfile(summaryDir, 'report.csv'));
disp(report(:, {'name', 'status', 'seconds'}));
fprintf('Reports: %s\n', summaryDir);
if any(report.status == "failed" | report.status == "blocked" | report.status == "completed_with_failures")
    warning('Thesis:IncompleteSuite', 'Some cases failed or are blocked; see report.csv and case logs.');
end
end
