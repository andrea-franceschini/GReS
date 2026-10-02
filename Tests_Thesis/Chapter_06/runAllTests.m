function report = runAllTests(varargin)
% Run Chapter 6 serially using the common runner.
suiteRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(suiteRoot, 'Utils'));
report = thesisRun(suiteRoot, 'Chapter', 6, varargin{:});
end
