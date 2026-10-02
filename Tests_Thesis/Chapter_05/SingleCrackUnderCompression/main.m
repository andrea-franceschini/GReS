function result = main(varargin)
% Test case 2: Single crack under compression (thesis section 5.4.3.2).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "SingleCrackUnderCompression", varargin{:});
end

function result = runExperiment(options)
% Graded native mesh; finite crack of half-length 1, inclination 20 degrees.
% Outside the crack the two independent blocks are tied with mortar.
a = cosd(20);
nFault = 60;
nOuter = 24;
nY = 45;
x = [linspace(-50, -a, nOuter + 1), linspace(-a, a, nFault + 1), linspace(a, 50, nOuter + 1)];
x = unique(x);
% Resolve the near-crack layer without filling the far field with tiny cells.
eta = linspace(0, 1, nY + 1).^2;
domains = Discretizer.empty;
for d = 1:2
    g = structuredMesh(x, eta, [-0.5 0.5]);
    c = g.coordinates;
    yf = c(:, 1) * tand(20);
    if d == 1
        % eta=0 is the crack, eta=1 the lower outer boundary.
        % Reverse logical y instead to maintain positive cell orientation.
        g = structuredMesh(x, 1 - fliplr(eta), [-0.5 0.5]);
        c = g.coordinates;
        yf = c(:, 1) * tand(20);
        g.coordinates(:, 2) = -50 + c(:, 2) .* (50 + yf);
        interfaceTag = 4;
    else
        g.coordinates(:, 2) = yf + c(:, 2) .* (50 - yf);
        interfaceTag = 3;
    end
    faces = find(g.surfaces.tag == interfaceTag);
    faceX = mean(reshape(g.coordinates(g.surfaces.connectivity(faces, :)', 1), 4, []), 1)';
    g.surfaces.tag(faces(abs(faceX) < a - 1e-10)) = 7;
    bc = Boundaries(g);
    zNodes = unique(g.surfaces.connectivity(ismember(g.surfaces.tag, [1 2]), :));
    thesisAddBC(bc, "planeStrain", "dirichlet", "node", "displacements", zNodes, 3, 0);
    thesisAddBC(bc, "leftLoad", "neumann", "surface", "displacements", find(g.surfaces.tag == 5), 1, 100);
    thesisAddBC(bc, "rightLoad", "neumann", "surface", "displacements", find(g.surfaces.tag == 6), 1, -100);
    if d == 1
        coords = g.coordinates;
        anchor = find(abs(coords(:, 1)) < 1e-10 & abs(coords(:, 2) + 50) < 1e-10);
        thesisAddBC(bc, "translation", "dirichlet", "node", "displacements", anchor, [1 2], 0);
        rotationAnchor = find(abs(coords(:, 1) + 50) < 1e-10 & abs(coords(:, 2) + 50) < 1e-10);
        thesisAddBC(bc, "rotation", "dirichlet", "node", "displacements", rotationAnchor, 2, 0);
    end
    domains(d) = Discretizer('grid', g, 'boundaries', bc, 'materials', Materials('Input/materials.xml'));
    domains(d).addPhysicsSolver('Poromechanics');
end
inp = struct('masterDomain', 2, 'slaveDomain', 1, 'masterSurface', 3, 'slaveSurface', 4, ...
           'multiplierType', "P0", 'Quadrature', struct('type', "SegmentBasedQuadrature", 'gaussOrder', 4));
interfaces = InterfaceSolver.add('MeshTying', domains, inp);
inp.masterSurface = 7;
inp.slaveSurface = 7;
inp.Coulomb = struct('cohesion', 0, 'frictionAngle', 30);
interfaces = InterfaceSolver.add('SolidMechanicsContact', domains, interfaces, inp);
out = OutState('outputFile', 'Output/singleCrack', 'matFileName', 'Output/history', 'printTimes', 1);
solver = NonLinearImplicit('simulationparameters', thesisStaticParameters(), 'domains', domains, 'interface', interfaces, 'output', out);
solver.simulationLoop();
contact = interfaces{2};
s = contact.getState();
c = contact.grids(MortarSide.slave).surfaces.center;
xi = c(:, 1) / cosd(20);
gt = reshape(s.tangentialGap, 2, [])';
gt = vecnorm(gt, 2, 2);
tn = s.traction(1:3:end);
K = 4 * (1 - 0.25^2) * 100 * sind(20) * (cosd(20) - sind(20) * tand(30)) / 15000;
exact = K * sqrt(max(0, 1 - xi.^2));
exactTn = -100 * sind(20)^2;
[~, perm] = sort(xi);
result.profile = table(xi(perm), gt(perm), exact(perm), tn(perm), ...
                     'VariableNames', {'xi', 'tangentialGap', 'analyticalGap', 'normalTraction'});
result.relativeGapError = norm(gt - exact) / norm(exact);
result.normalTractionError = norm(tn - exactTn) / sqrt(numel(tn));
assert(all(isfinite([gt; tn])), 'Nonfinite crack solution.');
writetable(result.profile, 'Output/profile.csv');
fig = figure;
tiledlayout(1, 2);
nexttile;
plot(xi(perm), gt(perm), 'o', xi(perm), exact(perm), '-');
xlabel('xi');
ylabel('Tangential gap');
legend('Mortar', 'Analytical');
nexttile;
plot(xi(perm), tn(perm), 'o');
yline(exactTn);
xlabel('xi');
ylabel('Normal traction');
exportgraphics(fig, 'Output/profile.pdf', 'ContentType', 'vector');
end
