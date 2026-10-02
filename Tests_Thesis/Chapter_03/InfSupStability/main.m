function result = main(varargin)
% Lagrange multiplier spaces and inf-sup stability (thesis section 3.2).
if nargin == 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'SuiteRoot')
    result = runExperiment(varargin{1});
    return
end
suiteRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(suiteRoot, 'Utils'));
result = thesisRun(suiteRoot, 'Cases', "InfSupStability", varargin{:});
end

function result = runExperiment(options)
% Migrated from Chapter_03/infsup/infSupPoisson.m.
result = struct();
%% infSup3D_refined_and_poisson.m
% Inf-sup estimate on a 3D two-patch mesh-tying problem using GReS.
%
% The estimate is computed from the generalized eigenvalue problem
%
%     B A^{-1} B' lambda = beta_h^2 Q lambda,
%
% where A is the Dirichlet-reduced bulk Poisson matrix, B is the
% Dirichlet-reduced mesh-tying constraint matrix, and Q is the chosen
% multiplier norm matrix.
%
% MeshTying stores the mortar blocks as
%
%     M*u_master + D*u_slave = 0.
%
% Therefore, after removing Dirichlet columns,
%
%     B = [M(:,keepMaster), D(:,keepSlave)].
%
% For P0 multipliers, Q_mass is the slave-interface face-area mass matrix.
% For standard/dual nodal multipliers, Q_mass is extracted from the slave
% side mortar block, i.e. Q_mass = -D(:,slaveInterfaceDofs).
%
% For the H^{-1/2}-scaled multiplier norm, the matrix used in the
% generalized eigenproblem is
%
%     Q = h_s * Q_mass,
%
% where h_s is the slave-side interface grid size.  Set multiplierNorm =
% "L2" if you instead want the pure L2 multiplier norm.

%% User parameters
% Choose the refinement on the two sides explicitly here.
% The two vectors must have the same length.
% NmList  = [2 4 8 16];      % master elements per coordinate direction
% NnmList = [4 8 16 32];     % slave/non-mortar elements per coordinate direction

% swap
NnmList  = [2 4 8 16];      % master elements per coordinate direction
NmList = [4 8 16 32];     % slave/non-mortar elements per coordinate direction

multTypes = ["standard", "P0", "P0stab"];     % any subset of ["standard", "P0", "dual"]
multiplierNorm = "HminusHalfScaled"; % "HminusHalfScaled" or "L2"

variable = "u";
masterSurface = 1;
slaveSurface  = 2;
fixMasterTags = setdiff(1:6, masterSurface);
fixSlaveTags  = setdiff(1:6, slaveSurface);

modifyBoundaryMultipliers = true;    % relevant only for nodal multipliers

% Poisson multiplier comparison. This also works on non-matching meshes.
% Choose the master/slave grid sizes independently.
runPoissonComparison = false;
poissonNm  = 4;     % master elements per coordinate direction
poissonNnm = 8;     % slave/non-mortar elements per coordinate direction
poissonMultTypes = ["standard", "P0"];

% Output names
if sum(NnmList) > sum(NmList)
    infSupPdfName  = "Output/infSup_patch3D_slaveFiner.pdf";
else
    infSupPdfName  = "Output/infSup_patch3D_slaveCoarser.pdf";
end

multPdfName  = "Output/mult.pdf";

%% Plot defaults
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultLineLineWidth', 1.4);
set(groot, 'defaultLineMarkerSize', 6);

%% Basic checks
if numel(NmList) ~= numel(NnmList)
    error("NmList and NnmList must have the same length.");
end
if any(NmList < 1) || any(NnmList < 1)
    error("All mesh sizes must be positive integers.");
end

% The h-scaled L2 multiplier norm is a discrete proxy, not the full H^{-1/2} norm.
%% Inf-sup refinement loop

if options.Smoke
    NmList = NmList(1);
    NnmList = NnmList(1);
end
nRefs = numel(NnmList);
nTypes = numel(multTypes);
beta = nan(nRefs, nTypes);
beta2 = nan(nRefs, nTypes);
hSlave = 1 ./ NnmList(:);
rankB = nan(nRefs, nTypes);
nMult = nan(nRefs, nTypes);

for it = 1:nTypes
    multType = multTypes(it);
    fprintf("\n==============================\n");
    fprintf("Multiplier type: %s\n", multType);
    fprintf("==============================\n");

    for ir = 1:nRefs
        Nm  = NmList(ir);
        Nnm = NnmList(ir);

        out = runInfSupCase(Nm, Nnm, multType, variable, masterSurface, slaveSurface, ...
                            fixMasterTags, fixSlaveTags, modifyBoundaryMultipliers, ...
                            multiplierNorm);

        beta(ir, it)  = out.beta;
        beta2(ir, it) = out.beta2;
        rankB(ir, it) = out.rankB;
        nMult(ir, it) = out.nMult;

        fprintf("Nm = %4d, Nnm = %4d, h_s = %.4e, rank(B) = %d/%d, beta = %.6e\n", ...
                Nm, Nnm, out.hSlave, out.rankB, out.nMult, out.beta);
    end
end

%% Plot inf-sup constant versus log(h_s)
fig = figure('Color', 'w');
hold on;
grid on;
box on;
markers = {'o', 's', '^', 'd', 'v', '>'};
for it = 1:nTypes
    mk = markers{1 + mod(it - 1, numel(markers))};
    plot(log(hSlave), beta(:, it), ['-' mk], 'DisplayName', sprintf('%s', multTypes(it)));
end
xlabel('$\log(h^{(1)})$');
ylabel('$\beta^*_h$');
legend('Location', 'best');
ylim([0 1]);
exportgraphics(fig, infSupPdfName, 'ContentType', 'vector');

%% Optional Poisson multiplier comparison on a possibly non-matching interface mesh
if runPoissonComparison
    if numel(poissonMultTypes) ~= 2 || ~all(ismember(["standard", "P0"], poissonMultTypes))
        error('poissonMultTypes must contain exactly "standard" and "P0" for the requested comparison.');
    end

    prof = struct();
    for it = 1:numel(poissonMultTypes)
        mt = poissonMultTypes(it);
        prof.(matlab.lang.makeValidName(mt)) = solvePoissonLambdaProfile( ...
            poissonNm, poissonNnm, mt, variable, masterSurface, slaveSurface, ...
            fixMasterTags, fixSlaveTags, modifyBoundaryMultipliers, multiplierNorm);
    end

    fig = figure('Color', 'w');
    hold on;
    grid on;
    box on;
    plot(prof.standard.s, prof.standard.lambda, 'k-o', 'DisplayName', 'standard');
    plot(prof.P0.s, prof.P0.lambda, 'r-s', 'DisplayName', 'P0');
    xlabel('$s$ along interface diagonal');
    ylabel('$\lambda$');
    legend('Location', 'best');
    exportgraphics(fig, multPdfName, 'ContentType', 'vector');
end

%% Store output
results.NmList = NmList(:);
results.NnmList = NnmList(:);
results.hSlave = hSlave;
results.multTypes = multTypes;
results.multiplierNorm = multiplierNorm;
results.beta = beta;
results.beta2 = beta2;
results.rankB = rankB;
results.nMult = nMult;
save('Output/infSup3D_results.mat', 'results');
result = results;

%% ========================================================================
%% Local functions
%% ========================================================================

end

function out = runInfSupCase(Nm, Nnm, multType, variable, masterSurface, slaveSurface, ...
                             fixMasterTags, fixSlaveTags, modifyBoundaryMultipliers, ...
                             multiplierNorm)
[A, B, Qmass, interface, domainMaster, domainSlave, keepMaster, keepSlave] = ...
    buildPatchMatrices(Nm, Nnm, multType, variable, masterSurface, slaveSurface, ...
                       fixMasterTags, fixSlaveTags, modifyBoundaryMultipliers);

hSlave = 1 / Nnm;
Q = applyMultiplierNormScaling(Qmass, hSlave, multiplierNorm);

% Schur complement in multiplier space.  Never form inv(Q)*S.
S = full(B * (A \ B'));
S = 0.5 * (S + S');

if strcmp(multType, "P0stab")
    H = full(interface.stabilizationMat);
    S = S + H;
    % H = interface.
end

M = full(Q \ S);
eigVals = eig(M);
eigVals = real(eigVals);

% Clean tiny roundoff below zero, but keep true zero modes.
tolEig = 1e-12 * max(1, max(abs(eigVals)));
eigVals(abs(eigVals) < tolEig) = 0;

beta2 = min(eigVals);
beta  = sqrt(max(beta2, 0));

out.beta = beta;
out.beta2 = beta2;
out.eigVals = eigVals;
out.hSlave = hSlave;
out.A = A;
out.B = B;
out.Qmass = Qmass;
out.Q = Q;
out.S = S;
out.rankB = rank(full(B));
out.nMult = size(B, 1);
out.interface = interface;
out.domainMaster = domainMaster;
out.domainSlave = domainSlave;
out.keepMaster = keepMaster;
out.keepSlave = keepSlave;
end

function profile = solvePoissonLambdaProfile(Nm, Nnm, multType, variable, masterSurface, slaveSurface, ...
                                             fixMasterTags, fixSlaveTags, modifyBoundaryMultipliers, ...
                                             multiplierNorm)
% The Poisson solve is the same saddle-point mesh-tying problem used for
% the inf-sup test and therefore works also for non-matching interface
% meshes. The extracted lambda profiles simply live on different
% abscissae for standard nodal and P0 face multipliers.

[A, B, ~, interface, domainMaster, domainSlave, keepMaster, keepSlave, T] = ...
    buildPatchMatrices(Nm, Nnm, multType, variable, masterSurface, slaveSurface, ...
                       fixMasterTags, fixSlaveTags, modifyBoundaryMultipliers);

nmult = size(B, 1);
K = [A, B'; B, sparse(nmult, nmult)];
rhs = [domainMaster.rhs{1}(keepMaster); ...
       domainSlave.rhs{1}(keepSlave); ...
       zeros(nmult, 1)];

sol = K \ rhs;
lambda = sol(end - nmult + 1:end);

profile = extractLambdaOnInterfaceDiagonal(interface, multType, lambda, T, Nnm);
end

function [A, B, Qmass, interface, domainMaster, domainSlave, keepMaster, keepSlave, T] = ...
          buildPatchMatrices(Nm, Nnm, multType, variable, masterSurface, slaveSurface, ...
                             fixMasterTags, fixSlaveTags, modifyBoundaryMultipliers)
% Put the master block above the slave block so that the selected
% interface surfaces are geometrically coincident.
gMaster = structuredMesh(Nm, Nm, Nm, [0 1], [0 1], [0.5 1]);
gSlave  = structuredMesh(Nnm, Nnm, Nnm, [0 1], [0 1], [0 0.5]);

% Manufactured field.  For -Delta u = f, the consistent forcing is
% f = 2*pi^2*y*sin(pi*x)*sin(pi*z).
uex = @(x, y, z) y .* sin(pi * x) .* sin(pi * z);
fex = @(x, y, z) -2 * pi^2 * y .* sin(pi * x) .* sin(pi * z);

bcMaster = Boundaries(gMaster);
bcMaster.addBC('name', "fixMaster", ...
               'type', "dirichlet", ...
               'field', "surface", ...
               'variable', variable, ...
               'entityListType', "tag", ...
               'entityList', fixMasterTags);
bcMaster.addBCEvent("fixMaster", 'time', 0.0, 'value', 0);

bcSlave = Boundaries(gSlave);
bcSlave.addBC('name', "fixSlave", ...
              'type', "dirichlet", ...
              'field', "surface", ...
              'variable', variable, ...
              'entityListType', "tag", ...
              'entityList', fixSlaveTags);
bcSlave.addBCEvent("fixSlave", 'time', 0.0, 'value', 0);

domainMaster = Discretizer('boundaries', bcMaster, 'grid', gMaster);
domainMaster.addPhysicsSolver('Poisson');
domainSlave = Discretizer('boundaries', bcSlave, 'grid', gSlave);
domainSlave.addPhysicsSolver('Poisson');

% Use analytical source if the local Poisson solver supports it.
try
    domainMaster.getPhysicsSolver("Poisson").setAnalSolution(uex, fex);
    domainSlave.getPhysicsSolver("Poisson").setAnalSolution(uex, fex);
catch
    warning('Could not set analytical Poisson solution/source. Continuing with default RHS.');
end

domainMaster.initialize();
domainSlave.initialize();
domains = [domainMaster; domainSlave];

if strcmp(multType, "P0stab")
    multType = "P0";
end

interfData = struct('masterDomain', 1, ...
                    'slaveDomain', 2, ...
                    'masterSurface', masterSurface, ...
                    'slaveSurface', slaveSurface, ...
                    'multiplierType', multType);
interfaces = InterfaceSolver.add("MeshTying", domains, interfData);
interface = interfaces{1};

% Initialize interface state consistently with the domains.
d = {domainMaster, domainSlave, interface};
for i = 1:numel(d)
    di = d{i};
    state = di.getState;
    di.setStateInit(state);
    di.setStateOld(state);
end

assembleDomain(domainMaster);
assembleDomain(domainSlave);
interface.assembleConstraint();

AfullMaster = domainMaster.J{1, 1};
AfullSlave  = domainSlave.J{1, 1};

dirMaster = bcMaster.getDofs("fixMaster", domainMaster.dofm);
dirSlave  = bcSlave.getDofs("fixSlave", domainSlave.dofm);

keepMaster = true(size(AfullMaster, 1), 1);
keepSlave  = true(size(AfullSlave, 1), 1);
keepMaster(dirMaster) = false;
keepSlave(dirSlave)   = false;

A = blkdiag(AfullMaster(keepMaster, keepMaster), ...
            AfullSlave(keepSlave, keepSlave));
A = sparse(0.5 * (A + A'));

M = interface.M;
D = interface.D;
Qmass = computeMultiplierMassMatrix(interface, multType);
T = speye(size(Qmass, 1));

if modifyBoundaryMultipliers && string(multType) ~= "P0"
    T = boundaryMultiplierTransformation(interface, bcSlave, "fixSlave");
    M = T * M;
    D = T * D;
    Qmass = T * Qmass * T';
end

B = [M(:, keepMaster), D(:, keepSlave)];
Qmass = sparse(0.5 * (Qmass + Qmass'));
end

function assembleDomain(domain)
% GReS versions differ on whether dt is required by assembleSystem.
try
    domain.assembleSystem();
catch
    try
        domain.assembleSystem(1.0);
    catch
        assembleSystem(domain, 1.0);
    end
end
end

function Q = applyMultiplierNormScaling(Qmass, hSlave, multiplierNorm)
switch string(multiplierNorm)
    case "HminusHalfScaled"
        Q = hSlave * Qmass;
    case "L2"
        Q = Qmass;
    otherwise
        error('Unknown multiplierNorm: %s. Use "HminusHalfScaled" or "L2".', multiplierNorm);
end
end

function Q = computeMultiplierMassMatrix(interf, multType)
s = MortarSide.slave;
nComp = interf.domains(s).dofm.getNumberOfComponents(interf.coupledVariables);

switch string(multType)
    case "P0"
        % P0 multipliers live on slave interface faces.  Their L2 mass
        % matrix is diagonal, with entries equal to the support area.
        area = interf.grids(s).surfaces.area(:);
        Q = spdiags(repelem(area, nComp), 0, nComp * numel(area), nComp * numel(area));

    case {"standard", "dual"}
        % For scalar standard multipliers, N_mult = N_slave interface
        % dofs, hence -D(:,slaveInterfaceDofs) is the consistent nodal
        % mass matrix.  For dual multipliers, MeshTying usually returns
        % the diagonal/biorthogonal version of this block.
        slaveDofs = getSlaveInterfaceDofs(interf);
        Q = -interf.D(:, slaveDofs);
        Q = sparse(0.5 * (Q + Q'));

    otherwise
        error("Unknown multiplier type: %s", multType);
end
end

function slaveDofs = getSlaveInterfaceDofs(interf)
s = MortarSide.slave;
gS = interf.grids(s);
dofmS = interf.domains(s).dofm;
fldS = dofmS.getVariableId(interf.coupledVariables);

nodes = unique(gS.surfaces.loc2glob(:), 'stable');
nodes = nodes(nodes > 0);
slaveDofs = dofmS.getLocalDoF(fldS, nodes);
end

function T = boundaryMultiplierTransformation(interf, bcSlave, bcName)
% Row transformation for nodal multipliers. Boundary multiplier rows are
% assigned to the nearest interior slave-interface node and then removed.
% This is the 3D analogue of endpoint row condensation in 2D mortar tests.

s = MortarSide.slave;
gS = interf.grids(s);
nComp = interf.domains(s).dofm.getNumberOfComponents(interf.coupledVariables);

nodes = unique(gS.surfaces.loc2glob(:), 'stable');
nodes = nodes(nodes > 0);
X = interf.domains(s).grid.coordinates(nodes, :);

bcNodes = bcSlave.getTargetEntities(bcName);
isBnd = ismember(nodes, bcNodes);

intNodeIds = find(~isBnd);
bndNodeIds = find(isBnd);

if isempty(bndNodeIds)
    T = speye(nComp * numel(nodes));
    return
end
if isempty(intNodeIds)
    error("All slave interface multiplier nodes are on the boundary.");
end

nOld = nComp * numel(nodes);
nNew = nComp * numel(intNodeIds);

rows = [];
cols = [];
vals = [];

% Keep all interior-node multiplier rows.
for k = 1:numel(intNodeIds)
    oldNode = intNodeIds(k);
    for c = 1:nComp
        rows(end + 1, 1) = nComp * (k - 1) + c; %#ok<AGROW>
        cols(end + 1, 1) = nComp * (oldNode - 1) + c; %#ok<AGROW>
        vals(end + 1, 1) = 1.0; %#ok<AGROW>
    end
end

% Add each removed boundary-node row to the closest interior node row.
Xi = X(intNodeIds, :);
for ib = reshape(bndNodeIds, 1, [])
    dx = Xi - X(ib, :);
    [~, nearLoc] = min(sum(dx.^2, 2));
    newNode = nearLoc;
    for c = 1:nComp
        rows(end + 1, 1) = nComp * (newNode - 1) + c; %#ok<AGROW>
        cols(end + 1, 1) = nComp * (ib - 1) + c; %#ok<AGROW>
        vals(end + 1, 1) = 1.0; %#ok<AGROW>
    end
end

T = sparse(rows, cols, vals, nNew, nOld);
end

function profile = extractLambdaOnInterfaceDiagonal(interface, multType, lambda, T, ~)
s = MortarSide.slave;
gs = interface.grids(s);

switch string(multType)
    case {"standard", "dual"}
        % Recover a full nodal vector before plotting, because boundary
        % multiplier rows may have been condensed away by T.

        mults = zeros(gs.nNodes, 1);
        innerMult = find(ismember(gs.surfaces.loc2glob, find(keepSlave)));

        mults(innerMult) = mult;

        Nid = find(abs(gsint.coordinates(:, 1) - gsint.coordinates(:, 2)) < 1e-6);

        profile.x = X(id, 1);
        profile.s = sqrt(2) * profile.x;
        profile.lambda = lambdaFull(loc);

    case "P0"
        C = gs.surfaces.center;
        tol = 1e-10;
        id = find(abs(C(:, 1) - C(:, 2)) < tol);
        [~, ord] = sort(C(id, 1));
        id = id(ord);
        profile.x = C(id, 1);
        profile.s = sqrt(2) * profile.x;
        profile.lambda = lambda(id);

    otherwise
        error("Unknown multiplier type: %s", multType);
end

% Defensive sorting and duplicate removal, useful for some grid orderings.
[profile.s, ord] = sort(profile.s(:));
profile.x = profile.x(ord);
profile.lambda = profile.lambda(ord);
end
