clc
clear
close all

% 1) homogeneous domain, single phase flow
% 2) time varying boundary conditions
% 3) heterogenous materials
% 4) initial (hydrostatic) condition
% 5) coupling with mechanics (one way/monolithic/fixed stress)

% create the mesh
grid = structuredMesh(11,11,7,[0 1000],[0 1000],[0 100]);

% a mesh with refinement, all using the internal grid generator
%grid = structuredMesh([4 7 4],[4 7 4],[3 5 3],[0 400 600 1000],[0 400 600 1000],[0 40 60 100]);
%grid = Grid();
%grid.importMesh('domain.vtk');
% show grid (no geometry processed)

grid.processGeometry();
% show grid (geometry processed)


% option 1: script based
materials = Materials();
%materials.addFluid('dynamicViscosity',1.1574e-11,'specificWeight',0.0,'compressibility',0.0);
materials.addFluid('dynamicViscosity',1.1574e-11,'specificWeight',0.0,'compressibility',4.4e-7);
materials.addSolid('name',"sand",'cellTags',1);
materials.addPorousRock("sand","permeability",1e-12,"porosity",0.375,'compressibility',1e-5);
%materials.addConstitutiveLaw("sand","Elastic",'youngModulus',1e3,'poissonRatio',0.25);

% option 2: xml input file
%materials2 = Materials('materials.xml');



% boundary conditions
bc = Boundaries(grid);

% pick the surfaces at the top (z = 10)
topSurfId = find(grid.surfaces.center(:,3) > 99.99);

bc.addBC('name',"topPressure",...
        'type',"dirichlet",...
        'field',"surface",...
        'entityListType',"bcList",...
        'entityList',topSurfId,...
        'variable',"pressure");

% pressure remains zero from 0 onward
bc.addBCEvent("topPressure",'time',0.0,'value',0.0);

% prescribe a source to a cell in the middle of the domain (extraction well)
% use the cell centroid to find the cell in the middle
wellCellsId = find(...
grid.cells.center(:,1) > 490.0 & grid.cells.center(:,1) < 510.0 & ...
grid.cells.center(:,2) > 490.0 & grid.cells.center(:,2) < 510.0 & ...  
grid.cells.center(:,3) > 48.0 & grid.cells.center(:,3) < 51.0);

bc.addBC('name',"extractionWell",...
        'type',"source",...
        'field',"cell",...
        'entityListType',"bcList",...
        'entityList',wellCellsId,...
        'variable',"pressure");
 
% this input is 1/TIME for a well prescribed on a cell!
vol = sum(grid.cells.volume(wellCellsId));
wellSrc = -1e4/vol;       % we prescribe an outflux of 1e3 m^3/year

bc.addBCEvent("extractionWell",'time',0.0,'value',wellSrc);

% discretizer
domain = Discretizer('Boundaries', bc, ...
                     'Materials',  materials, ...
                     'Grid',       grid);

domain.addPhysicsSolver("SinglePhaseFlowFVTPFA",'steadyState',0);

% simulation parameters
input = struct('Start',0.0,'End',10.0,'DtInit',2e-1,'DtMax',2e-1,'DtMin',2e-1,'incrementFactor',1.0);
simParam = SimulationParameters(input);


% input for the simulation output
out = OutState('printTimes',[0 0.1 0.25 0.5 1:10],'outputFile',"Output/flow",'matFileName',"Output/flow");

% set the solution scheme
solver = NonLinearImplicit('simulationparameters',simParam,'domains',domain, 'output', out);

% run the simulation
solver.simulationLoop();


%% Adding a less permeable material in the bottom of the domain and prescribing initial hydrostatic distribution

% set initial hydrostatic pressure
gamma = materials.getFluid.getSpecificWeight;

% set a cell tag to a bottom layer
bottomLayer = find(grid.cells.center(:,3) < 30.0);
grid.cells.tag(bottomLayer) = 2;
grid.cells.nTag = 2;

materials.addSolid('name',"clay",'cellTags',2);
materials.addPorousRock("clay","permeability",1e-15,"porosity",0.3,'compressibiliWHybrid TPFA
ty',1e-5);

domain = Discretizer('Boundaries', bc, ...
                     'Materials',  materials, ...
                     'Grid',       grid);
domain.addPhysicsSolver("SinglePhaseFlowFVTPFA",'steadyState',0);

domain.setState(gamma*(10-grid.cells.center(:,3)),"pressure");


out = OutState('printTimes',[0 0.1 0.25 0.5 1:3],'outputFile',"Output/flowClay",'matFileName',"Output/flow");


solver = NonLinearImplicit('simulationparameters',simParam,'domains',domain, 'output', out);

% run the simulation
solver.simulationLoop();

%% Mechanics (one way coupling)

% constitutive laws
materials.addConstitutiveLaw("sand","Elastic",'youngModulus',1e4,'poissonRatio',0.25);
materials.addConstitutiveLaw("clay","Elastic",'youngModulus',0.5e4,'poissonRatio',0.25);

% boundary condition for mechanics
bc.addBC('name',"fixBottom",...
        'type',"dirichlet",...
        'field',"surface",...
        'entityListType',"tags",...
        'entityList',1,...
        'variable',"displacements", ...
        'components',["x","y","z"]);

% pressure remains zero from 0 onward
bc.addBCEvent("fixBottom",'time',0.0,'value',0.0);


% pressure passed to mechanical model using dedicated boundary condition
bc.addBC('name',"pressureChange",...
        'type',"volumeForce",...
        'field',"cell",...
        'entityListType',"bcList",...
        'entityList',1:grid.cells.num,...
        'variable',"displacements");

flowModel = load('Output\flow.mat');
t = [flowModel.output.time];
p = [flowModel.output.pressure];

% impose a pressure change
deltaP = p - p(:,1);



% fix the bottom of the domain
for i = 1:numel(t)
  bc.addBCEvent("pressureChange",'time',t(i),'value',deltaP(:,i));
end

domain = Discretizer('Boundaries', bc, ...
                     'Materials',  materials, ...
                     'Grid',       grid);



domain.addPhysicsSolver("Poromechanics");

%domain.setState(gamma*(100-grid.cells.center(:,3)),"pressure");

out = OutState('printTimes',[0 0.1 0.25 0.5 1:3],'outputFile',"Output/coupledProblemOneWay",'matFileName',"Output/coupledProblem");

solver = NonLinearImplicit('simulationparameters',simParam,'domains',domain, 'output', out);

% run the simulation
solver.simulationLoop();



%% Adding mechanics


% pressure remains zero from 0 onward
bc.removeBC("pressureChange");

domain = Discretizer('Boundaries', bc, ...
                     'Materials',  materials, ...
                     'Grid',       grid);



domain.addPhysicsSolver("BiotFullyCoupled");

%domain.setState(gamma*(100-grid.cells.center(:,3)),"pressure");

out = OutState('printTimes',[0 0.1 0.25 0.5 1:3],'outputFile',"Output/coupledProblemFC",'matFileName',"Output/coupledProblem");

solver = NonLinearImplicit('simulationparameters',simParam,'domains',domain, 'output', out);

% run the simulation
solver.simulationLoop();


%% Fixed stress split

domain = Discretizer('Boundaries', bc, ...
                     'Materials',  materials, ...
                     'Grid',       grid);



domain.addPhysicsSolver("BiotFixedStressSplit");

%domain.setState(gamma*(100-grid.cells.center(:,3)),"pressure");

out = OutState('printTimes',[0 0.1 0.25 0.5 1:3],'outputFile',"Output/coupledProblemFS",'matFileName',"Output/coupledProblem");

solver = FixedStressSplit('simulationparameters',simParam,'domains',domain, 'output', out,'reltolerance',1e-2);

% run the simulation
solver.simulationLoop();













