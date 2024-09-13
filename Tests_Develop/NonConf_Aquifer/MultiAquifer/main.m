% MULTI-AQUIFER MODEL: Non conforming grid blocks

%%
clear
close all

%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('MULTIAQUIFER MODEL \n')
fprintf('___________________\n\n')
%% Define mesh objects
[d1,d2,d3,d4,d5] = deal(Mesh(),Mesh(),Mesh(),Mesh(),Mesh());
d1.importVTKmesh('model-2km/output_001.vtk');
d2.importVTKmesh('model-2km/output_002.vtk');
d3.importVTKmesh('model-2km/output_003.vtk');
d4.importVTKmesh('model-2km/output_004.vtk');
d5.importVTKmesh('model-2km/output_005.vtk');

%%
% n1 = find(all([d2.coordinates(:,1)==0, d2.coordinates(:,2)==0],2));
% n2 = find(all([d2.coordinates(:,1)==2000, d2.coordinates(:,2)==2000],2));
% n3 = find(all([d4.coordinates(:,1)==0, d4.coordinates(:,2)==0],2));
% n4 = find(all([d4.coordinates(:,1)==2000, d4.coordinates(:,2)==2000],2));

%% Utils
% bot_nodes = load('bot_nodes');
% writeBCfiles('press_L2','NodeBC','Dir','Flow',[],'injection_L2',[0 1 5 6 100],[0 100 100 0 0],[403 405]);
% writeBCfiles('press_L4','NodeBC','Dir','Flow',[],'withdraw_L4',[0 1 5 6 100],[0 -100 -100 0 0],[402 404]);
% writeBCfiles('bot_fix_L5','NodeBC','Dir','Poro',["x","y","z"],'fixed_bottom',0,0,bot_nodes);
%% Setup multidomain simulation
simParam = SimulationParameters('simParam.dat');
model = buildModelStruct('domains2km.dat',simParam);

%% Setup computations between non conforming domains
% add manually the interfaces between domains
% model(1).Grid.topology.addSurface(1,'Mesh/vol_1_surf_2.topol');
% model(2).Grid.topology.addSurface(1,'Mesh/vol_2_surf_1.topol');
% model(2).Grid.topology.addSurface(2,'Mesh/vol_2_surf_2.topol');
% model(3).Grid.topology.addSurface(1,'Mesh/vol_3_surf_1.topol');
% model(3).Grid.topology.addSurface(2,'Mesh/vol_3_surf_2.topol');
% model(4).Grid.topology.addSurface(1,'Mesh/vol_4_surf_1.topol');
% model(4).Grid.topology.addSurface(2,'Mesh/vol_4_surf_2.topol');
% model(5).Grid.topology.addSurface(1,'Mesh/vol_5_surf_1.topol');
model(1).Grid.topology.addSurface(1,'model-2km/vol_1_surf_2.topol');
model(2).Grid.topology.addSurface(1,'model-2km/vol_2_surf_1.topol');
model(2).Grid.topology.addSurface(2,'model-2km/vol_2_surf_2.topol');
model(3).Grid.topology.addSurface(1,'model-2km/vol_3_surf_1.topol');
model(3).Grid.topology.addSurface(2,'model-2km/vol_3_surf_2.topol');
model(4).Grid.topology.addSurface(1,'model-2km/vol_4_surf_1.topol');
model(4).Grid.topology.addSurface(2,'model-2km/vol_4_surf_2.topol');
model(5).Grid.topology.addSurface(1,'model-2km/vol_5_surf_1.topol');

%% Load mesh glue structure (avoid recomputing mortar every time)
fprintf('Computing mortar utils \n')
mG = MeshGlue(model,'interfaces.dat');
fprintf('Done computing mortar utils\n')

%meshGlueStruct = load("meshGlue_aquiferSlave.mat");
%mG = meshGlueStruct.mG;
%clear meshGlueStruct


%% modify interpolation operator to be sparse
for i = 1:numel(mG.interfaces)
    E = mG.interfaces(i).InterpOperator;
    E(abs(E) < 1e-4) = 0;
    E = E./sum(E,2);
    mG.interfaces(i).InterpOperator = sparse(E);
end

fprintf('Model Setup completed \n')
% fName = 'interfaces_slaveFlow.dat';
% mG = MeshGlue(model,fName);
% fprintf('Done reading interfaces \n')
%% MULTIDOMAIN SOLUTION STRATEGY
NSolv_MD = NonLinearSolverMultiDomain(simParam,model,mG);
% Solution loop
NSolv_MD.NonLinearLoop_multiAquifer();

%%
for i = 1:numel(model) 
   mG.model(i).OutState.finalize()
end


