% MULTI-AQUIFER MODEL: Non conforming grid blocks

%%
clear
close all

%% Define mesh objects
[d1,d2,d3,d4,d5] = deal(Mesh(),Mesh(),Mesh(),Mesh(),Mesh());
d1.importVTKmesh('Mesh/output_001.vtk');
d2.importVTKmesh('Mesh/output_002.vtk');
d3.importVTKmesh('Mesh/output_003.vtk');
d4.importVTKmesh('Mesh/output_004.vtk');
d5.importVTKmesh('Mesh/output_005.vtk');
%% Utils
% bot_nodes = load('bottom_nodes.dat');
% writeBCfiles('press_inj_L2','NodeBC','Dir','Flow',[],'injection_L2',[0 1 5 6 100],[0 100 100 0 0],[603 606]);
% writeBCfiles('press_erog_L4','NodeBC','Dir','Flow',[],'withdraw_L4',[0 1 5 6 100],[0 -100 -100 0 0],[601 603]);
% writeBCfiles('fix','NodeBC','Dir','Poro',["x","y","z"],'fixed_bottom',0,0,bot_nodes);
%% Setup multidomain simulation
simParam = SimulationParameters('simParam.dat');
model = buildModelStruct('domains.dat',simParam);

%% Setup computations between non conforming domains
% add manually the interfaces between domains
model(1).Grid.topology.addSurface(1,'Mesh/vol_1_surf_2.topol');
model(2).Grid.topology.addSurface(1,'Mesh/vol_2_surf_1.topol');
model(2).Grid.topology.addSurface(2,'Mesh/vol_2_surf_2.topol');
model(3).Grid.topology.addSurface(1,'Mesh/vol_3_surf_1.topol');
model(3).Grid.topology.addSurface(2,'Mesh/vol_3_surf_2.topol');
model(4).Grid.topology.addSurface(1,'Mesh/vol_4_surf_1.topol');
model(4).Grid.topology.addSurface(2,'Mesh/vol_4_surf_2.topol');
model(5).Grid.topology.addSurface(1,'Mesh/vol_5_surf_1.topol');

%% Load mesh glue structure (avoid recomputing mortar every time)
meshGlueStruct = load("meshGlue_aquiferSlave.mat");
mG = meshGlueStruct.mG;
clear meshGlueStruct

% fName = 'interfaces_slaveFlow.dat';
% mG = MeshGlue(model,fName);
% fprintf('Done reading interfaces \n')
%% MULTIDOMAIN SOLUTION STRATEGY
NSolv_MD = NonLinearSolverMultiDomain(simParam,model,mG);
% Solution loop
NSolv_MD.NonLinearLoop();

for i = 1:numel(model)
%% 
   model(i).OutState.finalize()
end


