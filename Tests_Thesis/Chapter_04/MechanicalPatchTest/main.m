function result = main(varargin)
% Mechanical patch test (supplementary) (thesis section 4.4.3).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "MechanicalPatchTest", varargin{:});
end

function result = runExperiment(options)
% Supplementary compression patch formerly in Chapter_04/PatchTest3D.
% Zero Poisson ratio matches a uniaxial constant-strain patch with fixed base.
methods = ["SegmentBasedQuadrature", "ElementBasedQuadrature", "RBFquadrature"];
result = struct();
result.methods = methods;
result.maxDisplacementError = zeros(3, 1);
for k = 1:numel(methods)
    domains = Discretizer.empty;
    for d = 1:2
        n = 4;
        if d == 2
            n = 3;
        end
        g = structuredMesh(n, n, n, [0 1], [0 1], [d - 1 d]);
        bc = Boundaries(g);
        if d == 1
            nodes = unique(g.surfaces.connectivity(g.surfaces.tag == 1, :));
            thesisAddBC(bc, "base", "dirichlet", "node", "displacements", nodes, [1 2 3], 0);
        else
            thesisAddBC(bc, "load", "neumann", "surface", "displacements", find(g.surfaces.tag == 2), 3, -0.5);
        end
        mat = Materials();
        mat.addSolid('name', "solid", 'cellTags', 1);
        mat.addConstitutiveLaw("solid", "Elastic", 'youngModulus', 1000, 'poissonRatio', 0);
        domains(d) = Discretizer('grid', g, 'materials', mat, 'boundaries', bc);
        domains(d).addPhysicsSolver('Poromechanics');
    end
    inp = struct('masterDomain', 1, 'slaveDomain', 2, 'masterSurface', 2, 'slaveSurface', 1, ...
     'multiplierType', "P0", 'Quadrature', struct('type', methods(k), 'gaussOrder', 5, 'nInt', 4));
    interfaces = InterfaceSolver.add('MeshTying', domains, inp);
    out = OutState('outputFile', fullfile('Output', methods(k)), 'printTimes', 1);
    solver = NonLinearImplicit('simulationparameters', thesisStaticParameters(), 'domains', domains, 'interface', interfaces, 'output', out);
    solver.simulationLoop();
    err = 0;
    for d = 1:2
        disp = reshape(domains(d).getState('displacements'), 3, [])';
        exact = zeros(size(disp));
        exact(:, 3) = -0.5 / 1000 * domains(d).grid.coordinates(:, 3);
        err = max(err, max(abs(disp - exact), [], 'all'));
    end
    result.maxDisplacementError(k) = err;
    if k == 1
        assert(err < 1e-8, 'Thesis:MechanicalPatch', 'Mechanical patch failed.');
    end
end
end
