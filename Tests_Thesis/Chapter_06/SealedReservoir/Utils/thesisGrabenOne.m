function result = thesisGrabenOne(options)
% Migrated from Chapter_06/Graben/grabenThesis.m.
result = struct();
% this model represents a field with multiple fracture
% two main fracture bound an aquifer and are explicitely represented as the
% domain of two mortar discretizations
% smaller fractures are simulated using EFEM
% in the aquifer, Biot model is used for mechanical/solid coupling

%% generate grid and skew external surfaces

input = readInput('Input/input_SS.xml');

% g1 = structuredMesh(15,10,[4 2 4],[0 2e3],[0 3e3],[0 450 550 1e3]);
% %g2 = structuredMesh(10,10,[6 10 6],[2e3 3e3],[0 3e3],[0 450 550 1e3]);
% g2 = structuredMesh(5,[5 5 5],[10 6 10],[2e3 3e3],[0 1e3 2e3 3e3],[0 450 550 1e3]);
% g3 = structuredMesh(15,10,[4 2 4],[3e3 5e3],[0 3e3],[0 450 550 1e3]);

% old converged graben_fine_01
g1 = structuredMesh(4, 6, 6, [0 1.3e3], [0 3e3], [0 1e3]);
g2 = structuredMesh(8, [10 12 10], [10 12 10], [1.3e3 2e3], [0 1e3 2e3 3e3], [0 400 600 1e3]);
g3 = structuredMesh(10, [12 14 12], [12 14 12], [2e3 3e3], [0 1e3 2e3 3e3], [0 430 630 920]);
g4 = structuredMesh(8, [10 12 10], [10 12 10], [3e3 3.7e3], [0 1e3 2e3 3e3], [0 400 600 1e3]);
g5 = structuredMesh(4, 6, 6, [3.7e3 5e3], [0 3e3], [0 1e3]);

grids = [g1; g2; g3; g4; g5];

% Number of embedded fractures in grids 1 through 5.  The strict geometry
% rules below require zero fractures in grids 1 and 5 and at most 30 total.
fracturesPerDomain = [0 8 10 8 0];
[domains, fracs] = setDomains(grids, fracturesPerDomain);

% Same struct array produced by readInput for repeated <Fracture ... />
% entries; it can be passed directly to the embedded-fracture solver.
solBase = input.Solver;
solBase.BiotFullyCoupled.Mechanics.Fracture = fracs;

simparam = SimulationParameters(input.SimulationParameters);

doms = [1 2; 2 3; 4 3; 5 4];
surfs = [6 5; 6 5; 5 6; 5 6];

% stress initialization
for i = 1:numel(domains)
    domains(i).addPhysicsSolver("PreStressMechanics");
end

outPre = OutState('outputFile', 'preStress', 'printTimes', 0.025, 'solvePrintTimes', 0);

interfaces = {};

for i = 1:size(doms, 1)
    interfInput = struct('masterDomain', doms(i, 1), ...
                         'slaveDomain', doms(i, 2), ...
                         'masterSurface', surfs(i, 1), ...
                         'slaveSurface', surfs(i, 2), ...
                         'multiplierType', "P0", ...
                         'stabilizationScale', 1.0);

    interfaces = InterfaceSolver.add("MeshTying", domains, interfaces, interfInput);

end

solver = NonLinearImplicit('simulationparameters', simparam, 'domains', domains, 'output', outPre, 'interface', interfaces);

solver.simulationLoop();

stress = cell(numel(domains), 1);
for i = 1:numel(stress)
    stress{i} = getState(domains(i), "stress");
end

%% coupled model with pressure

c = g3.cells.center;
isAquifer = c(:, 3) > 400 & c(:, 3) < 600 & c(:, 2) > 1e3 & c(:, 2) < 2e3;
g3.cells.tag(isAquifer) = 2;
g3.cells.nTag = g3.cells.nTag + 1;

% domains = setDomains(g1,g2,g3);

well = c(:, 3) > 500 & c(:, 3) < 520 & ...
  c(:, 1) > 2.4e3 & c(:, 1) < 2.6e3 & ...
  c(:, 2) > 1.4e3 & c(:, 2) < 1.6e3;
wellId = find(well);
assert(~isempty(wellId), 'Thesis:NoWellCells', 'No cells selected for the well.');

domains(3).bcs.addBC('name', "wellPress", ...
                     'type', "dirichlet", ...
                     'field', "cell", ...
                     'variable', "pressure", ...
                     'entityListType', "bclist", ...
                     'entityList', wellId);
domains(3).bcs.addBCEvent("wellPress", 'time', 0.0, 'value', 0.0);
% domains(2).bcs.addBCEvent("wellPress",'time',0.4,'value',-3e3);
% domains(2).bcs.addBCEvent("wellPress",'time',0.405,'value',0.0);
domains(3).bcs.addBCEvent("wellPress", 'time', 1.0, 'value', -3.5e3);
domains(3).bcs.addBCEvent("wellPress", 'time', 2.0, 'value', 0.0);
% domains(2).bcs.addBCEvent("wellPress",'time',1.0,'value',-2.5e3);
% domains(2).bcs.addBCEvent("wellPress",'time',2.0,'value',-100);
% domains(2).bcs.addBCEvent("wellPress",'time',10.0,'value',-100.0);

domains(3).materials.addSolid('name', "sandAquifer", "specificWeight", 21.0, 'cellTags', 2);
domains(3).materials.addPorousRock("sandAquifer", 'permeability', 1e-13);
domains(3).materials.addConstitutiveLaw("sandAquifer", "Elastic", 'youngModulus', 0.85e6, 'poissonRatio', 0.25);

for i = 1:numel(domains)
    domains(i) = Discretizer('Boundaries', domains(i).bcs, ...
                             'Materials', domains(i).materials, ...
                             'Grid', grids(i));
end

domains(1).addPhysicsSolver("Poromechanics");
if options.WithSubsidiaryFractures
    n = 0;
    for i = 2:4
        sol = solBase;
        sol.BiotFullyCoupled.Mechanics.Fracture = fracs(n + 1:n + fracturesPerDomain(i));
        if i == 3
            domains(i).addPhysicsSolver("BiotFullyCoupled", sol.BiotFullyCoupled);
        else
            domains(i).addPhysicsSolver("EFEMaugmented", sol.BiotFullyCoupled.Mechanics);
        end
        n = n + fracturesPerDomain(i);
    end
else
    domains(2).addPhysicsSolver("Poromechanics");
    domains(3).addPhysicsSolver("BiotFullyCoupled", struct("Flow", struct('targetRegions', 2)));
    domains(4).addPhysicsSolver("Poromechanics");
end
domains(5).addPhysicsSolver("Poromechanics");

for i = 1:numel(domains)
    domains(i).setState(stress{i}, "stress");
end

outSim = OutState('outputFile', 'grabenHorst', 'printTimes', 0:0.1:2.0, 'solvePrintTimes', 1);

input.SimulationParameters.End = 2.0;
simparam = SimulationParameters(input.SimulationParameters);

% plotMesh(g1,'test1');
% plotMesh(g2,'test2');
% plotMesh(g3,'test3');

interfaces = InterfaceSolver.addInterfaces(domains, input.Interface);

solver = NonLinearImplicit('simulationparameters', simparam, 'domains', domains, 'output', outSim, 'interface', interfaces);

solver.simulationLoop();

result.fractures = fracs;
result.withSubsidiaryFractures = options.WithSubsidiaryFractures;
result.wellCells = wellId;

%% materials, boundary conditions

end

function [domains, fracs] = setDomains(grids, fracturesPerDomain)

% Return domains without any physics solver applied, together with a flat
% struct array matching the repeated Fracture field read from the XML.
%
% fracturesPerDomain can be either:
%   [nGrid1 nGrid2 nGrid3 nGrid4 nGrid5], or
%   [nGrid2 nGrid3 nGrid4].
% Grids 1 and 5 are outside the prescribed damage zones and must contain no
% embedded fractures.  Grid 3 fractures are divided between its two faults.

if numel(grids) ~= 5
    error("setDomains:GridCount", ...
          "The graben model requires the five grids g1 through g5.");
end

if ~isnumeric(fracturesPerDomain) || ~isreal(fracturesPerDomain) || ...
    ~isvector(fracturesPerDomain)
    error("setDomains:InvalidFractureCounts", ...
          "fracturesPerDomain must be a real numeric vector.");
end

fracturesPerDomain = fracturesPerDomain(:).';
if numel(fracturesPerDomain) == 3
    fracturesPerDomain = [0 fracturesPerDomain 0];
elseif numel(fracturesPerDomain) ~= numel(grids)
    error("setDomains:InvalidFractureCounts", ...
          "Use either [nGrid2 nGrid3 nGrid4] or one count for each of the five grids.");
end

if any(~isfinite(fracturesPerDomain)) || ...
    any(fracturesPerDomain < 0) || ...
    any(fracturesPerDomain ~= fix(fracturesPerDomain))
    error("setDomains:InvalidFractureCounts", ...
          "Every fracture count must be a finite non-negative integer.");
end
if any(fracturesPerDomain([1 5]) ~= 0)
    error("setDomains:FractureLocation", ...
          "Embedded fractures are allowed only in grids 2, 3, and 4.");
end
if sum(fracturesPerDomain) > 30
    error("setDomains:TooManyFractures", ...
          "The total number of embedded fractures must not exceed 30.");
end

% Points on the unskewed main-fault faces.  Keeping these before skewGrid is
% useful because each requested normal defines the final planar face through
% the centre of the original boundary.
pFault2 = faceCentre(grids(2), 6);
pFault3West = faceCentre(grids(3), 5);
pFault3East = faceCentre(grids(3), 6);
pFault4 = faceCentre(grids(4), 5);

skew = struct([]);
skew(1).surfTag = 6;
skew(1).normal  = [1 0.1 0.3];
skewGrid(grids(2), skew);

skew = struct([]);
skew(1).surfTag = 5;
skew(1).normal  = [-1 -0.1 -0.3];
skew(2).surfTag = 6;
skew(2).normal  = [1 -0.05 -0.2];
skewGrid(grids(3), skew);

skew = struct([]);
skew(1).surfTag = 5;
skew(1).normal  = [-1 0.05 0.2];
skewGrid(grids(4), skew);

for i = 1:numel(grids)
    processGeometry(grids(i));
end

% The requested grid-3 population is split as evenly as possible between its
% two faults.  Generation uses MATLAB's current random stream; call rng(seed)
% before setDomains when an exactly reproducible realization is required.
fracs = makeDamageZoneFractures(grids, ...
                                pFault2, pFault3West, pFault3East, pFault4, fracturesPerDomain);

% tag aquifer

mat = Materials();
mat.addSolid('name', "rock", "specificWeight", 21.0, 'cellTags', 1);
mat.addConstitutiveLaw("rock", "Elastic", 'youngModulus', 1e6, 'poissonRatio', 0.25, 'bulkModulusRatio', 0.25);
mat.addFluid('compressibility', 4.59e-7, 'dynamicViscosity', 3.1608e-14);

matWeak = Materials();
matWeak.addSolid('name', "rock", "specificWeight", 21.0, 'cellTags', 1);
matWeak.addConstitutiveLaw("rock", "Elastic", 'youngModulus', 1e6, 'poissonRatio', 0.25);
matWeak.addFluid('compressibility', 4.59e-7, 'dynamicViscosity', 3.1608e-14);

% old material set
% mat = Materials();
% mat.addSolid('name',"rock","specificWeight",21.0,'cellTags',1);
% mat.addConstitutiveLaw("rock","Elastic",'youngModulus',1.0e6,'poissonRatio',0.25);
% mat.addFluid('compressibility',4.59e-7,'dynamicViscosity',3.1608e-14);
%
%
% matWeak = Materials();
% matWeak.addSolid('name',"rock","specificWeight",21.0,'cellTags',1);
% matWeak.addConstitutiveLaw("rock","Elastic",'youngModulus',1e6,'poissonRatio',0.25);
% matWeak.addFluid('compressibility',4.59e-7,'dynamicViscosity',3.1608e-14);
% matWeak.addSolid('name',"sandAquifer","specificWeight",21.0,'cellTags',2);
% matWeak.addConstitutiveLaw("sandAquifer","Elastic",'youngModulus',1e5,'poissonRatio',0.25);

domainList = cell(numel(grids), 1);
for i = 1:numel(grids)
    bound = Boundaries(grids(i));
    bound.addBC('name', "z_fix", ...
                'type', "dirichlet", ...
                'field', "surface", ...
                'variable', "displacements", ...
                'entityListType', "tag", ...
                'entityList', 1, ...
                'components', "z");
    bound.addBCEvent("z_fix", 'time', 0.0, 'value', 0.0);

    bound.addBC('name', "y_fix", ...
                'type', "dirichlet", ...
                'field', "surface", ...
                'variable', "displacements", ...
                'entityListType', "tag", ...
                'entityList', [3 4], ...
                'components', "y");
    bound.addBCEvent("y_fix", 'time', 0.0, 'value', 0.0);

    % Only the two external x faces are fixed; the internal faces participate
    % in the two mortar interfaces.
    if i == 1 || i == numel(grids)
        xTag = 5 + (i == numel(grids));
        bound.addBC('name', "x_fix", ...
                    'type', "dirichlet", ...
                    'field', "surface", ...
                    'variable', "displacements", ...
                    'entityListType', "tag", ...
                    'entityList', xTag, ...
                    'components', "x");
        bound.addBCEvent("x_fix", 'time', 0.0, 'value', 0.0);
    end

    mat = Materials();
    mat.addSolid('name', "rock", "specificWeight", 21.0, 'cellTags', 1);
    mat.addConstitutiveLaw("rock", "Elastic", 'youngModulus', 1e6, 'poissonRatio', 0.25, 'bulkModulusRatio', 0.25);
    mat.addFluid('compressibility', 4.59e-7, 'dynamicViscosity', 3.1608e-14);

    material = mat;

    if i == 3
        material = matWeak;
    end

    domainList{i} = Discretizer('Boundaries', bound, ...
                                'Materials', material, ...
                                'Grid', grids(i));
end

domains = vertcat(domainList{:});

end

function fracs = makeDamageZoneFractures(grids, ...
                                         pFault2, pFault3West, pFault3East, pFault4, fracturesPerDomain)
% MAKEDAMAGEZONEFRACTURES Generate complete rectangular subsidiary faults.
%
% The fracture origin is assumed to be the rectangle centre.  Candidates
% are accepted only when all four corners lie inside the host-grid bounding
% box and on the host-grid side of the corresponding main-fault plane.

nGrid3West = ceil(fracturesPerDomain(3) / 2);
nGrid3East = floor(fracturesPerDomain(3) / 2);

% Outward normals of the four main-fault faces.  The two descriptions of
% each main fault remain separate because the neighbouring grids may be
% nonconforming.
populations = struct( ...
                     'gridId', {2, 3, 3, 4}, ...
                     'count', {fracturesPerDomain(2), nGrid3West, nGrid3East, ...
           fracturesPerDomain(4)}, ...
                     'point', {pFault2, pFault3West, pFault3East, pFault4}, ...
                     'normal', {[1 0.1 0.3], [-1 -0.1 -0.3], ...
            [1 -0.05 -0.2], [-1 0.05 0.2]});

emptyFrac = struct('origin', [], 'normal', [], 'lengthVec', [], ...
                   'widthVec', [], 'dimensions', [], 'cohesion', 0.0, 'frictionAngle', 28.0);
fracs = repmat(emptyFrac, sum([populations.count]), 1);

maxAttempts = 500;
k = 0;

for p = 1:numel(populations)
    grid = grids(populations(p).gridId);
    xyzMin = min(grid.coordinates, [], 1);
    xyzMax = max(grid.coordinates, [], 1);
    domainSize = xyzMax - xyzMin;

    baseNormal = populations(p).normal;
    baseNormal = baseNormal / norm(baseNormal);
    faultPoint = populations(p).point;

    % The subsidiary faults in the reference sketch are shorter than the
    % main fault, steep, and approximately subparallel to it.  Their centres
    % are concentrated near the fault core but fill the damage-zone width.
    zoneThickness = min(450, 0.55 * domainSize(1));
    lengthRange = [0.14 0.24] * domainSize(2);  % strike direction
    widthRange  = [0.24 0.42] * domainSize(3);  % down-dip direction

    % Keep the full rectangle slightly away from both the outer box and the
    % main-fault plane.  This avoids apparent cuts caused by roundoff or by a
    % rectangle lying exactly on a domain boundary.
    clearance = max(1, 0.01 * min(domainSize));

    for j = 1:populations(p).count
        k = k + 1;
        accepted = false;

        for attempt = 1:maxAttempts
            % Small deviations preserve the subparallel subsidiary-fault pattern
            % while avoiding an artificial set of perfectly coplanar rectangles.
            [normal, lengthVec, widthVec] = ...
              perturbedFaultFrame(baseNormal, 8, 10);

            dimensions = [randomInRange(lengthRange), ...
                          randomInRange(widthRange)];

            % Draw a point on the main-fault plane.  Sampling away from the y/z
            % ends improves acceptance without determining validity by itself;
            % validity is established below from the four actual corners.
            y = xyzMin(2) + (0.12 + 0.76 * rand) * domainSize(2);
            z = xyzMin(3) + (0.12 + 0.76 * rand) * domainSize(3);
            x = faultPoint(1) - ...
              (baseNormal(2) * (y - faultPoint(2)) + ...
                         baseNormal(3) * (z - faultPoint(3))) / baseNormal(1);
            pointOnFault = [x y z];

            % The rectangle is tilted relative to the main fault.  Consequently,
            % its centre must be at least this far inside the domain to keep every
            % corner on the correct side of the main-fault plane.
            projectedHalfSize = 0.5 * ( ...
                                     dimensions(1) * abs(dot(lengthVec, baseNormal)) + ...
                                     dimensions(2) * abs(dot(widthVec, baseNormal)));
            minDepth = projectedHalfSize + clearance;

            if minDepth >= zoneThickness
                continue
            end

            if j == 1
                % Closest complete rectangle: nearly tangent to the main fault, but
                % never centred on it (which would place half the rectangle outside).
                depth = minDepth;
            else
                % rand^1.7 biases the population toward the fault core, as expected
                % for subsidiary faults in a damage zone.
                depth = minDepth + (zoneThickness - minDepth) * rand^1.7;
            end

            origin = pointOnFault - depth * baseNormal;
            corners = rectangleCorners(origin, lengthVec, widthVec, dimensions);

            insideBox = all(corners >= xyzMin + clearance, 'all') && ...
              all(corners <= xyzMax - clearance, 'all');
            insideFault = all((corners - faultPoint) * baseNormal' <= -clearance);

            if insideBox && insideFault
                accepted = true;
                break
            end
        end

        if ~accepted
            error('makeDamageZoneFractures:NoAdmissibleRectangle', ...
                  ['Could not place a complete fracture in grid %d after %d ', ...
                             'attempts. Reduce the fracture-size ranges or the clearance.'], ...
                  populations(p).gridId, maxAttempts);
        end

        fracs(k).origin = origin;
        fracs(k).normal = normal;
        fracs(k).lengthVec = lengthVec;
        fracs(k).widthVec = widthVec;
        fracs(k).dimensions = dimensions;
        fracs(k).cohesion = 0.0;
        fracs(k).frictionAngle = 28.0;
    end
end

end

function [normal, lengthVec, widthVec] = perturbedFaultFrame(baseNormal, ...
                                                           maxNormalDeviation, maxStrikeDeviation)

% Reference strike and down-dip directions of the main fault.
strike = [0 1 0] - dot([0 1 0], baseNormal) * baseNormal;
strike = strike / norm(strike);
dip = cross(baseNormal, strike);
dip = dip / norm(dip);

% Uniform azimuth and area-uniform sampling in a narrow cone around the
% main-fault normal produce steep, subparallel subsidiary faults.
azimuth = 2 * pi * rand;
deviation = deg2rad(maxNormalDeviation) * sqrt(rand);
normal = cos(deviation) * baseNormal + ...
  sin(deviation) * (cos(azimuth) * strike + sin(azimuth) * dip);
normal = normal / norm(normal);

% Reproject the reference strike into the perturbed plane and apply a small
% in-plane rotation, yielding elongated but non-identical rectangles.
strike = strike - dot(strike, normal) * normal;
strike = strike / norm(strike);
width0 = cross(normal, strike);
width0 = width0 / norm(width0);
strikeRotation = deg2rad(maxStrikeDeviation) * (2 * rand - 1);
lengthVec = cos(strikeRotation) * strike + sin(strikeRotation) * width0;
lengthVec = lengthVec / norm(lengthVec);
widthVec = cross(normal, lengthVec);
widthVec = widthVec / norm(widthVec);

end

function corners = rectangleCorners(origin, lengthVec, widthVec, dimensions)
halfLength = 0.5 * dimensions(1) * lengthVec;
halfWidth = 0.5 * dimensions(2) * widthVec;
corners = [origin + halfLength + halfWidth; ...
           origin + halfLength - halfWidth; ...
           origin - halfLength + halfWidth; ...
           origin - halfLength - halfWidth];
end

function value = randomInRange(range)
value = range(1) + diff(range) * rand;
end

function centre = faceCentre(grid, tag)
xyzMin = min(grid.coordinates, [], 1);
xyzMax = max(grid.coordinates, [], 1);
centre = 0.5 * (xyzMin + xyzMax);
if tag == 5
    centre(1) = xyzMin(1);
elseif tag == 6
    centre(1) = xyzMax(1);
else
    error("faceCentre:UnsupportedTag", "Only x-face tags 5 and 6 are used.");
end
end
