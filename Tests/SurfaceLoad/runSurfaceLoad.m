clear
close all
clc


sim = struct('Start',0.0,...
  'End',10,...
  'DtInit',0.25,...
  'DtMax',0.25,...
  'DtMin',1e-5,...
  'RelativeTolerance',1.e-8,...
  'AbsoluteTolerance',1.e-9, ...
  'LinearSolver',struct('useMatlab',0));


simParam = SimulationParameters(sim);

grid = Grid();
grid.importMesh("Test01_tetra.msh")
grid.processGeometry;

%
%----------------------------- MATERIALS -----------------------------
%

% materials
mat = Materials();
mat.addFluid('dynamicViscosity',1.157e-11,'specificWeight',0.0,'compressibility',4.59e-7);

names = ["Silt","LowerSiltyClay","Sand","UpperSiltyClay"];
for i = 1:4
  mat.addSolid('name',names(i),'cellTags',i);
end

mat.addPorousRock("Silt","permeability",1e-13,"porosity",0.3);
mat.addConstitutiveLaw("Silt","Elastic",'youngModulus',2.129e3,'poissonRatio',0.3);

mat.addPorousRock("LowerSiltyClay","permeability",1e-14,"porosity",0.3);
mat.addConstitutiveLaw("LowerSiltyClay","Elastic",'youngModulus',2.129e3,'poissonRatio',0.3);

mat.addPorousRock("Sand","permeability",[1e-11,1e-11,1e-11],"porosity",0.3);
mat.addConstitutiveLaw("Sand","Elastic",'youngModulus',4.503e3,'poissonRatio',0.3);

mat.addPorousRock("UpperSiltyClay","permeability",1e-14,"porosity",0.3);
mat.addConstitutiveLaw("UpperSiltyClay","Elastic",'youngModulus',2.129e3,'poissonRatio',0.3);

%
%------------------------ BOUNDARY CONDITIONS ------------------------
%
bc = Boundaries(grid);

% dirichlet displacements fixed
bc.addBC('name',"dir_Poro_fixed",...
        'type',"dirichlet",...
        'field',"node",...
        'entityListType',"bclistfile",...
        'entityList',"dir_Poro_fixed/list",...
        'variable',"displacements");

bc.addBCEvent("dir_Poro_fixed",'time',0.0,'value',0.0);

% surface top load
bc.addBC('name',"neuSurf_Poro_topLoad",...
        'type',"neumann",...
        'field',"surface",...
        'entityListType',"bclistfile",...
        'entityList',"neuSurf_Poro_TopLoad/list",...
        'variable',"displacements", ...
        'components',"z");

bc.addBCEvent("neuSurf_Poro_topLoad",'time',0.0,'value',0.0);
bc.addBCEvent("neuSurf_Poro_topLoad",'time',3.0,'value',-8.0);
bc.addBCEvent("neuSurf_Poro_topLoad",'time',100.0,'value',-8.0);

% symmetry bcs
bc.addBC('name',"dir_Poro_Sym",...
        'type',"dirichlet",...
        'field',"node",...
        'entityListType',"bclistfile",...
        'entityList',"dir_Poro_Sym/list",...
        'variable',"displacements");

bc.addBCEvent("dir_Poro_Sym",'time',0.0,'value',0.0);


% top surface drained (zero pressure)
bc.addBC('name',"dir_Flow_TopDrained",...
        'type',"dirichlet",...
        'field',"node",...
        'entityListType',"bclistfile",...
        'entityList',"dir_Flow_TopDrained/list",...
        'variable',"pressure");

bc.addBCEvent("dir_Flow_TopDrained",'time',0.0,'value',0.0);

% Set the input file
%fileName = ["dir_Poro_fixed.dat","neuSurf_Poro_topLoad.dat","dir_Poro_Sym.dat","neuSurf_Flow_Impermeable.dat","dir_Flow_TopDrained.dat"];
%
%
% Create and set the print utility
out = OutState('outputFile',"surfaceLoad",'printTimes',0:0.5:10);

% discretizer
domain = Discretizer('Boundaries', bc, ...
                     'Materials',  mat, ...
                     'Grid',       grid);

domain.addPhysicsSolver("BiotFullyCoupled",'flowSolver',"SinglePhaseFlowFEM");


solver = NonLinearImplicit('simulationparameters',simParam,'domains',domain, 'output', out);
gresLog().setVerbosity(2);
solver.simulationLoop();


%%
% -------------------------- POST PROCESSING ------------------------------

nodes = [143;215;287;359];

disp = [out.results.displacements];
press = [out.results.pressure];
times = [out.results.time];

pressplot = press(nodes,:);
dispPlot = disp(nodes*3,:);
figure(1)
plot(times,pressplot,'k-o')
xlabel('Time (days)')
ylabel('Pressure (kPa)')
xlim([0 12])
ylim([-1 6])
grid on

figure(2)
plot(times,1e3*dispPlot,'k-o')
xlabel('Tempo (giorni)')
ylabel('Spostamenti verticali (mm)')
xlim([0 12])
ylim([-40 5])
grid on


timesInd = find(ismember([out.results.time],[1,2,3,5,10]));
nodes_subsidence = intersect(find(grid.coordinates(:,3) ==50),find(grid.coordinates(:,2) ==0));
x = sort(grid.coordinates(nodes_subsidence,1));
subs = sort(disp(3*nodes_subsidence,timesInd),'ascend');
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.')
figure(3)
plot(x,1000*subs);
xlabel('r (m)')
ylabel('Vertical displacements (mm)')
xlim([0 50])
ylim([-60 5])
time_string = "Day  " + times(timesInd);
legend(time_string)
grid on









