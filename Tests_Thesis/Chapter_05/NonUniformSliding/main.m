function result = main(varargin)
% Test case 1: Non uniform sliding (thesis section 5.5.4.1).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "NonUniformSliding", varargin{:});
end

function result = runExperiment(options)
% Migrated from Chapter_05/ConstantSlidingEFEM_borja/testConstantSlidingEFEM.m.
result = struct();

% Extract the directory containing the script

fname = 'Input/constantSlidingEFEM.xml';

params = readInput(fname);

simparams = SimulationParameters(params.SimulationParameters);

grid = structuredMesh(4, 2, 8, [0 2], [0, 0.5], [0 4]);

mat = Materials(params.Materials);

bc = Boundaries(grid, params.BoundaryConditions);
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('boundaries', bc, ...
                     'materials', mat, ...
                     'grid', grid);

domain.addPhysicsSolvers(params.Solver);

out = OutState('outputFile', "Output/EFEM_borja", 'printTimes', 1.0);

solver = NonLinearImplicit('simulationparameters', simparams, ...
                           'domains', domain, ...
                           'output', out);
solver.simulationLoop();
result.fractureJump = domain.getState("fractureJump");
assert(all(isfinite(result.fractureJump)));

% get tangential gap
% gt = abs(getState(domain,"fractureJump"));
% gt = gt(2:3:end);
% anGt = 0.1*sqrt(2);
% tol = 1e-6;
% assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")

%

end
