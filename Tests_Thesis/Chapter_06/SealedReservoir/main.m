function result = main(varargin)
% Sealed reservoir in a fractured graben-horst formation (thesis section 6.3.4).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "SealedReservoir", varargin{:});
end

function result = runExperiment(options)
base = pwd;
restore = onCleanup(@() cd(base)); %#ok<NASGU>
result = struct();
for withSubs = [true false]
    name = "withSubsidiaryFractures";
    if ~withSubs
        name = "withoutSubsidiaryFractures";
    end
    out = fullfile(base, 'Output', name);
    mkdir(out);
    copyfile(fullfile(base, 'Input'), fullfile(out, 'Input'));
    cd(out);
    options.WithSubsidiaryFractures = withSubs;
    result.(name) = thesisGrabenOne(options);
end
end
