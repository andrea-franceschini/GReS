function result = main(varargin)
% Standalone EFEM and mortar comparison, run serially after initGReS.
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
p = inputParser;
addParameter(p, 'Visible', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Seed', 1, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});
options = p.Results;
caseRoot = fileparts(mfilename('fullpath'));
gresRoot = getappdata(0, 'gres_root');
assert((ischar(gresRoot) || isstring(gresRoot)) && ...
       isfile(fullfile(gresRoot, 'initGReS.m')), ...
       'Thesis:GReSNotInitialized', 'Run initGReS once before the comparison.');
options.GReSRoot = gresRoot;
oldDir = pwd;
oldPath = path;
oldVisible = get(groot, 'defaultFigureVisible');
oldStream = RandStream.getGlobalStream;
oldStreamState = oldStream.State;
cleanup = onCleanup(@() restoreSession(oldDir, oldPath, oldVisible, oldStream, oldStreamState)); %#ok<NASGU>
addpath(fullfile(caseRoot, 'Utils'));
out = fullfile(caseRoot, 'Output');
if isfolder(out)
    rmdir(out, 's');
end
mkdir(out);
cd(caseRoot);
if ~options.Visible
    set(groot, 'defaultFigureVisible', 'off');
end
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', options.Seed));
diary(fullfile(out, 'run.log'));
result = runExperiment(options);
save(fullfile(out, 'result.mat'), 'result', 'options', '-v7.3');
fprintf('Comparison outputs: %s\n', out);
end

function restoreSession(folder, oldPath, visible, oldStream, oldStreamState)
diary off;
cd(folder);
path(oldPath);
oldStream.State = oldStreamState;
RandStream.setGlobalStream(oldStream);
set(groot, 'defaultFigureVisible', visible);
end

function result = runExperiment(options)
% Both formulations remain together, with independent input/cache/output trees.
base = pwd;
restore = onCleanup(@() cd(base)); %#ok<NASGU>
result = struct();
result.failedConfigurations = 0;
for model = ["EFEM", "Mortar"]
    modelDir = fullfile(base, 'Output', model);
    mkdir(modelDir);
    copyfile(fullfile(base, 'Input', model), fullfile(modelDir, 'Input'));
    cd(modelDir);
    data = feval("thesisCompare" + model, options);
    result.(model) = data;
    result.failedConfigurations = result.failedConfigurations + data.failedConfigurations;
end
end
