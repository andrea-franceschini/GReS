function result = main(varargin)
% Test case 1: Constant sliding (thesis section 5.4.3.1).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "ConstantSliding", varargin{:});
end

function result = runExperiment(options)
% Two blocks separated by y=0.7+x. Native mesh, no VTK/Gmsh dependency.
n = 12;
if options.Smoke
    n = 4;
end
domains = Discretizer.empty;
for d = 1:2
    g = structuredMesh(n, n, 2, [0 2], [0 1], [0 0.5]);
    c = g.coordinates;
    interfaceY = 0.7 + c(:, 1);
    if d == 1
        g.coordinates(:, 2) = c(:, 2) .* interfaceY;
    else
        g.coordinates(:, 2) = interfaceY + c(:, 2) .* (4 - interfaceY);
    end
    bc = Boundaries(g);
    planeNodes = unique(g.surfaces.connectivity(ismember(g.surfaces.tag, [1 2]), :));
    thesisAddBC(bc, "planeStrain", "dirichlet", "node", "displacements", planeNodes, 3, 0);
    if d == 1
        bottom = unique(g.surfaces.connectivity(g.surfaces.tag == 3, :));
        thesisAddBC(bc, "bottomY", "dirichlet", "node", "displacements", bottom, 2, 0);
        anchor = bottom(abs(g.coordinates(bottom, 1)) < 1e-10);
        thesisAddBC(bc, "anchorX", "dirichlet", "node", "displacements", anchor, 1, 0);
    else
        top = unique(g.surfaces.connectivity(g.surfaces.tag == 4, :));
        thesisAddBC(bc, "topY", "dirichlet", "node", "displacements", top, 2, -0.1);
    end
    domains(d) = Discretizer('grid', g, 'boundaries', bc, 'materials', Materials('Input/materials.xml'));
    domains(d).addPhysicsSolver('Poromechanics');
end
inp = struct('masterDomain', 1, 'slaveDomain', 2, 'masterSurface', 4, 'slaveSurface', 3, ...
 'multiplierType', "P0", 'Coulomb', struct('cohesion', 0, 'frictionAngle', 5.71), ...
 'Quadrature', struct('type', "SegmentBasedQuadrature", 'gaussOrder', 4));
interfaces = InterfaceSolver.add('SolidMechanicsContact', domains, inp);
out = OutState('outputFile', 'Output/constantSliding', 'matFileName', 'Output/history', 'printTimes', 1);
solver = NonLinearImplicit('simulationparameters', thesisStaticParameters(), 'domains', domains, 'interface', interfaces, 'output', out);
solver.simulationLoop();
s = interfaces{1}.getState();
gap = reshape(s.tangentialGap, 2, [])';
result.tangentialGap = sqrt(sum(gap.^2, 2));
result.expectedGap = 0.1 * sqrt(2);
result.maxError = max(abs(result.tangentialGap - result.expectedGap));
assert(all(isfinite(result.tangentialGap)), 'Nonfinite sliding solution.');
assert(result.maxError < 1e-6, 'Thesis:ConstantSliding', 'Constant-sliding analytical check failed.');
writetable(table(result.tangentialGap), 'Output/tangentialGap.csv');
end
