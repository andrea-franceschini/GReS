function result = main(varargin)
% Influence of the relative refinement (thesis section 3.3.1).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "InfluenceOfRelativeRefinement", varargin{:});
end

function result = runExperiment(options)
base = pwd;
restore = onCleanup(@() cd(base)); %#ok<NASGU>
result = struct();
for factor = [0.5 2]
    name = sprintf('slaveFactor_%g', factor);
    out = fullfile(base, 'Output', name);
    mkdir(out);
    copyfile(fullfile(base, 'Input'), fullfile(out, 'Input'));
    cd(out);
    options.SlaveFactor = factor;
    result.(matlab.lang.makeValidName(name)) = thesisRelativeRefinementOne(options);
end
end
