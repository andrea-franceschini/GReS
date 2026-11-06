clear
clc
close all
% White et al., CMAME 357 (2019), 112575: Section 6.1, Table 1, Fig. 1.
% Adapted from gresSPE10_mech.m. Run initGReS and start MRST beforehand.
%
% REPRODUCTION LIMITS
% - The three channel boxes below reconstruct Fig. 1; their internal face
%   coordinates are not tabulated in the paper. They are an explicit,
%   adjustable interpretation, not an author-supplied geometry dataset.
% - Like the supplied reference, this script uses uniform initial pressure
%   and zero gravity: an incremental, gravity-free benchmark variant. The
%   paper starts from hydrostatic/lithostatic equilibrium. Do not interpret
%   this script as reproducing the paper's absolute initial stress field.
% - MultiphasePoromechanics and PeacemanWells must be available from the
%   same GReS branch as the uploaded reference (not the public checkout).

mrstModule add ad-core ad-blackoil ad-props mrst-gui
gravity reset off
assert(exist('MultiphasePoromechanics','class') == 8, ...
    'Use the GReS branch containing MultiphasePoromechanics.');
assert(exist('PeacemanWells','class') == 8, ...
    'Use the GReS branch containing PeacemanWells.');

%% Mesh and reconstructed staircase channel
% Domain dimensions read from Fig. 1: 500 x 500 x 250 m.
% 32 x 32 x 16 cells give 88,307 unconstrained unknowns (3 nodal
% displacements + 2 cell variables), matching the paper's ~88k base mesh.
refinement = 1;                    % 1, 2, 4, ... for uniform refinement
assert(refinement >= 1 && refinement == fix(refinement));
n = refinement*[16,16,8];
boxSize = [500,500,250]*meter;
grid = structuredMesh(n(1),n(2),n(3), ...
    [0 boxSize(1)],[0 boxSize(2)],[0 boxSize(3)]);
grid.processGeometry();
xyz = grid.cells.center;

% Each row is [xmin xmax ymin ymax zmin zmax] in metres.
% Successive boxes share a vertical face over a nonzero height, so the
% channel is face-connected and does not depend on edge-only connections.
% All breaks align with the base mesh; edit here to refine the reconstruction.

% permeable regions
c = grid.cells.center;
xc = c(:,1);
yc = c(:,2);
zc = c(:,3);
Z = boxSize(3);

stair = false(grid.cells.num,4);
stair(:,1) = (xc > 0.0 & xc < 500.0) & (yc > 0.0 & yc < 250.0) & (zc > 0.0 & zc < 0.25*Z);
stair(:,2) = (xc > 0.0 & xc < 250.0) & (yc > 0.0 & yc < 500.0) & (zc > 0.25*Z & zc < 0.5*Z);
stair(:,3) = (xc > 0.0 & xc < 500.0) & (yc > 250.0 & yc < 500.0) & (zc > 0.5*Z & zc < 0.75*Z);
stair(:,4) = (xc > 250.0 & xc < 500.0) & (yc > 0.0 & yc < 500.0) & (zc > 0.75*Z & zc < Z);
isChannel = any(stair,2);

grid.cells.tag(:) = 2;              % low-permeability seal
grid.cells.tag(isChannel) = 1;      % permeable channel
assert(any(isChannel) && any(~isChannel));

%% Materials: Table 1, SI units
materials = Materials();
materials.addFluid();
regionNames = ["channel","seal"];
porosity = [0.20,0.05];
permeability = [1000,1]*milli*darcy; % isotropic, contrast 1000:1
for k = 1:2
    materials.addSolid('name',regionNames(k),'cellTags',k);
    materials.addPorousRock(regionNames(k), ...
        'permeability',permeability(k), ...
        'porosity',porosity(k), ...
        'biotCoefficient',1.0, ...
        'compressibility',0.0);     % incompressible grains, b = 1
    materials.addConstitutiveLaw(regionNames(k),"Elastic", ...
        'youngModulus',5000*mega*Pascal,'poissonRatio',0.25);
end
% No SPE10 data, burden layers, or extra empirical rock compressibility.
% With gravity off, the tabulated grain density (2650 kg/m^3) is inactive.

%% Two immiscible fluids: equations (41)--(42)
srW = 0.2; srNW = 0.2;
pRef = 6000*psia;                 % datum retained from supplied reference;
                                  % not a specified staircase datum in paper
muW = 0.3*centi*poise;
muNW = 3.0*centi*poise;
cW = 4.34e-4/(mega*Pascal);
cNW = 1.98e-4/(mega*Pascal);
rhoW0 = 1035*kilogram/meter^3;
rhoNW0 = 863*kilogram/meter^3;
rhoW = @(p) rhoW0*(1+cW*(p-pRef));
rhoNW = @(p) rhoNW0*(1+cNW*(p-pRef));
krW = coreyPhaseRelpermAD(2,srW,1,srW+srNW);
krNW = coreyPhaseRelpermAD(2,srNW,1,srW+srNW);
fluid.wet = struct('viscosity',muW,'compressibility',cW, ...
    'relPerm',krW,'density',rhoW);
fluid.nwet = struct('viscosity',muNW,'compressibility',cNW, ...
    'relPerm',krNW,'density',rhoNW);
fluid.rock = struct('refPressure',pRef);

%% Boundary conditions: impervious exterior, roller sides/bottom, free top
bc = Boundaries(grid);
% structuredMesh tags: 1 bottom, 2 top, 3/4 y-normal, 5/6 x-normal.
surfaceTags = [1,3,4,5,6];
components = [3,2,2,1,1];
for k = 1:numel(surfaceTags)
    name = "roller_"+string(surfaceTags(k));
    bc.addBC('name',name,'type',"dirichlet",'field',"surface", ...
        'entityListType',"tag",'entityList',surfaceTags(k), ...
        'variable',"displacements",'components',components(k));
    bc.addBCEvent(name,'time',0.0,'value',0.0);
end
% No flow boundary conditions => natural zero flux. No top displacement BC.

domain = Discretizer('Boundaries',bc,'Materials',materials,'Grid',grid);
% Both seal and channel participate in two-phase flow and mechanics.
domain.addPhysicsSolver("MultiphasePoromechanics", ...
    'Flow',struct('targetRegions',[1,2]));
twoph = getFlowSolver(domain.getPhysicsSolver("MultiphasePoromechanics"));
twoph.setFluidProperties(fluid);
domain.setState(pRef*ones(grid.cells.num,1),"pressure");
domain.setState(srW*ones(grid.cells.num,1),"saturation");

%% Wells: opposite ends, completed only through the channel
% x-fastest indexing, as in the uploaded reference. Select using coordinates
% to avoid any dependence on the ordering of cells within the mesh.
dx = boxSize(1)/n(1); dy = boxSize(2)/n(2);
prodId = find(c(:,1) > boxSize(1) - dx & c(:,2) < dy & c(:,3) < 0.25 * Z);
injId = find(c(:,1) > boxSize(1) - dx & c(:,2) < dy & c(:,3) > 0.75 * Z);
assert(~isempty(injId) && ~isempty(prodId),'Empty channel well completion.');
[~,order] = sort(xyz(injId,3)); injId = injId(order);
[~,order] = sort(xyz(prodId,3)); prodId = prodId(order);

rampTime = 1*day;
deltaP = 5*mega*Pascal;
% Linear ramp over day 1, then constant target BHP through day 100.
wellTimes = [0,rampTime,100*day];
wellMod = PeacemanWells(domain);
wellMod.addWell(injId,'control',"bhp", ...
    'time',wellTimes,'value',pRef+[0,deltaP,deltaP], ...
    'radius',0.1524*meter, ...
    'name','injector','dir','z','comp',[1,0],'sign',+1);
wellMod.addWell(prodId,'control',"bhp", ...
    'time',wellTimes,'value',pRef-[0,deltaP,deltaP], ...
    'radius',0.1524*meter, ...
    'name','producer','dir','z','sign',-1);
% Skin = 0 is assumed to be the default, as in the supplied reference API.
twoph.addWellModel(wellMod);

%% Time integration and output
% Initial/max/end times and Newton tolerance from Table 1. The exponential
% growth factor is not given in the paper; retain 1.1 from the reference.
simpar = struct('End',5*day,'DtInit',0.1*day,'DtMax',1*day, ...
    'DtMin',1e-5*day,'choppingFactor',2,'incrementFactor',1.1, ...
    'RelativeTolerance',1e-5);
simparams = SimulationParameters(simpar);
% The previous standalone CPR class supports flow-only systems. This model
% needs a three-field mechanics/flow preconditioner for iterative solves.
% Keep the reference script's solver selection here; do not force flow CPR.
if ~isfolder('Output'), mkdir('Output'); end
printUtils = OutState('printTimes',unique([0:day:simpar.End,92*day]), ...
    'outputFile',"Output/staircase_mech",'vtkFormat',"binary");

fprintf('Staircase: %d cells; %d channel cells; %d seal cells.\n', ...
    grid.cells.num,nnz(isChannel),nnz(~isChannel));
fprintf('Channel completions: injector %d, producer %d.\n', ...
    numel(injId),numel(prodId));
fprintf('Unconstrained DoFs: %d.\n',3*prod(n+1)+2*prod(n));
solver = NonLinearImplicit('domains',domain, ...
    'simulationparameters',simparams,'output',printUtils);
solver.simulationLoop();
