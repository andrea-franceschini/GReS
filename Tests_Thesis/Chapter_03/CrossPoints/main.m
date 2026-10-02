function result = main(varargin)
% Cross points (thesis section 3.3.3).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "CrossPoints", varargin{:});
end

function result = runExperiment(options)
% Migrated from Chapter_03/infsup/crossPointPatch.m.
result = struct();

usol = @(x, y, z) z;
fsol = @(x, y, z) 0 * x .* y .* z;

g = Grid.empty;

multType = "P0";
stabScale = 1.0;

nf = 3;
nc = 4;
% eight blocks
g(1) = structuredMesh(nf, nf, nf, [0, 0.5], [0 0.5], [0 0.5]);
g(2) = structuredMesh(nc, nc, nc, [0.5, 1], [0 0.5], [0 0.5]);
g(3) = structuredMesh(nf, nf, nf, [0.5, 1], [0.5 1], [0 0.5]);
g(4) = structuredMesh(nc, nc, nc, [0, 0.5], [0.5, 1], [0 0.5]);

g(5) = structuredMesh(nc, nc, nc, [0, 0.5], [0 0.5], [0.5, 1]);
g(6) = structuredMesh(nf, nf, nf, [0.5, 1], [0 0.5], [0.5, 1]);
g(7) = structuredMesh(nc, nc, nc, [0.5, 1], [0.5 1], [0.5, 1]);
g(8) = structuredMesh(nf, nf, nf, [0, 0.5], [0.5, 1], [0.5, 1]);

bc = Boundaries.empty;

for i = 1:8
    bc(i) = Boundaries(g(i));
end

bc(1) = Boundaries(g(1));
bc(1).addBC('name', "fix", ...
  'type', "dirichlet", ...
  'field', "node", ...
  'variable', "u", ...
  'entityListType', "bcList", ...
  'entityList', 1);
bc(1).addBCEvent("fix", 'time', 0.0, 'value', 0.0);

% neumann faces for each domain
neumannFaces = [3, 5
                6, 3
                4, 6
                4, 5];

for i = 1:4

    k = i;

    bc(i).addBC('name', "vertflux", ...
      'type', "neumann", ...
      'field', "surface", ...
      'variable', "u", ...
      'entityListType', "tags", ...
      'entityList', 1);

    bc(i).addBCEvent("vertflux", 'time', 0.0, 'value', -1);

    bc(i).addBC('name', "nf1", ...
      'type', "neumann", ...
      'field', "surface", ...
      'variable', "u", ...
      'entityListType', "tags", ...
      'entityList', neumannFaces(k, 1));

    bc(i).addBCEvent("nf1", 'time', 0.0, 'value', 0);

    bc(i).addBC('name', "nf2", ...
      'type', "neumann", ...
      'field', "surface", ...
      'variable', "u", ...
      'entityListType', "tags", ...
      'entityList', neumannFaces(k, 2));

    bc(i).addBCEvent("nf2", 'time', 0.0, 'value', 0);

end

k = 0;

for i = 5:8

    k = k + 1;

    bc(i).addBC('name', "vertflux", ...
      'type', "neumann", ...
      'field', "surface", ...
      'variable', "u", ...
      'entityListType', "tags", ...
      'entityList', 2);

    bc(i).addBCEvent("vertflux", 'time', 0.0, 'value', 1);

    bc(i).addBC('name', "nf1", ...
      'type', "neumann", ...
      'field', "surface", ...
      'variable', "u", ...
      'entityListType', "tags", ...
      'entityList', neumannFaces(k, 1));

    bc(i).addBCEvent("nf1", 'time', 0.0, 'value', 0);

    bc(i).addBC('name', "nf2", ...
      'type', "neumann", ...
      'field', "surface", ...
      'variable', "u", ...
      'entityListType', "tags", ...
      'entityList', neumannFaces(k, 2));

    bc(i).addBCEvent("nf2", 'time', 0.0, 'value', 0);
end

domains = Discretizer.empty;
for i = 1:8
    domains(i) = Discretizer('boundaries', bc(i), 'grid', g(i));
    domains(i).addPhysicsSolver('Poisson', 'gaussOrder', 4);
    domains(i).getPhysicsSolver("Poisson").setAnalSolution(usol, fsol);
end

% interfaces (reduntant pattern)
mortars = [ % bottom
           1, 2
           1, 4
           3, 2
           3, 4
           % top
           6, 5
           6, 7
           8, 5
           8, 7
           % bot2top
           1, 5
           2, 6
           3, 7
           4, 8
           ];

surfs =   [ % bottom
           6, 5
           4, 3
           3, 4
           5, 6
           % top
           5, 6
           4, 3
           3, 4
           6, 5
           % bot2top
           2, 1
           2, 1
           2, 1
           2, 1
           ];

% Interfaces forming a spanning tree:
% all eight domains are connected, with no redundant constraint cycle.

% 8 interfaces, still too many
% mortars = [
%   1,2;
%   1,4;
%   3,2;
%   3,4;
%   1,5;
%   2,6;
%   3,7;
%   4,8
%   ];
%
% surfs = [
%   6,5;   % 1 - 2
%   4,3;   % 1 - 4
%   3,4;   % 3 - 2
%   5,6;
%   2,1;   % 1 - 5
%   2,1;   % 6 - 2
%   2,1;   % 3 - 7
%   2,1    % 8 - 4
%   ];

% mortars = [
%     1,2;
%     1,4;
%     3,2;
%     1,5;
%     6,2;
%     3,7;
%     8,4
% ];
%
% surfs = [
%     6,5;   % 1 - 2
%     4,3;   % 1 - 4
%     3,4;   % 3 - 2
%     2,1;   % 1 - 5
%     1,2;   % 6 - 2
%     2,1;   % 3 - 7
%     1,2    % 8 - 4
% ];

interfaces = {};

for i = 1:size(mortars, 1)
    interfInput = struct('masterDomain', mortars(i, 1), ...
      'slaveDomain', mortars(i, 2), ...
      'masterSurface', surfs(i, 1), ...
      'slaveSurface', surfs(i, 2), ...
      'multiplierType', multType, ...
      'stabilizationScale', stabScale);

    interfaces = InterfaceSolver.add("MeshTying", domains, interfaces, interfInput);

end

outName = fullfile('Output', strcat("interf_", multType));
solver = NonLinearImplicit('simulationparameters', SimulationParameters('Input/simParam.xml'), ...
  'output', OutState('outputFile', outName, 'printTimes', 1), ...
  'domains', domains, ...
  'interface', interfaces);

solver.simulationLoop();
result.maxNodalError = 0;
for k = 1:numel(domains)
    result.maxNodalError = max(result.maxNodalError, max(abs(domains(k).getState('u') - domains(k).grid.coordinates(:, 3))));
end
result.multipliers = cellfun(@(x) x.getState('multipliers'), interfaces, 'UniformOutput', false);
assert(result.maxNodalError < 1e-8, 'Thesis:CrossPointPatch', 'Cross-point affine patch failed.');

end
