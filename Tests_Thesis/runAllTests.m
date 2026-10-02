function report = runAllTests(varargin)
% Run all thesis cases serially after initGReS.
% report = runAllTests('Smoke',true);
v = gresLog().getVerbosity();
gresLog().setVerbosity(-1);
suiteRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(suiteRoot, 'Utils'));
report = thesisRun(suiteRoot, varargin{:});
gresLog().setVerbosity(v);
end
