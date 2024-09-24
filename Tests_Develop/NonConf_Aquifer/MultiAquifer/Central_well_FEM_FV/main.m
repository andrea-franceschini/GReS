% MULTI-AQUIFER MODEL: Non conforming grid blocks
%%
clear
close all

%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('MULTIAQUIFER MODEL \n')
fprintf('___________________\n\n')
%% Define mesh objects
[d1,d2,d3,d4,d5] = deal(Mesh(),Mesh(),Mesh(),Mesh(),Mesh());
d1.importMesh('model-mixed/output_001.vtk');
d2.importMesh('model-mixed/output_002.vtk');
d3.importMesh('model-mixed/output_003.vtk');
d4.importMesh('model-mixed/output_004.vtk');
d5.importMesh('model-mixed/output_005.vtk');

% trans = 1000;
% C = {d1 d2 d3 d4 d5};
% totNod = 0;
% totCells = 0;
% for i = 1:numel(C)
% %     C{i}.coordinates(:,1) = C{i}.coordinates(:,1) - trans;
% %     plotFunction(C{i},strcat('dom_trans_',num2str(i)),zeros(C{i}.nNodes,1));
%     totNod = totNod+C{i}.nNodes;
%     totCells = totCells+C{i}.nCells;
% end

%save('mesh.mat',"d1","d2","d3","d4","d5");
%%
% n1 = find(all([d2.coordinates(:,1)==0, d2.coordinates(:,2)==0],2));
% n2 = find(all([d2.coordinates(:,1)==2000, d2.coordinates(:,2)==2000],2));
% n3 = find(all([d4.coordinates(:,1)==0, d4.coordinates(:,2)==0],2));
% n4 = find(all([d4.coordinates(:,1)==2000, d4.coordinates(:,2)==2000],2));
%bot_nodes = find(d5.coordinates(:,3)==-300);
% % 
%[~,pos] = min(sum([abs(d2.cellCentroid(:,1)), abs(d2.cellCentroid(:,2))],2));

%% Utils
% bot_nodes = load('bot_nodes');
% cL2inj = 9409:9412;
% cL2erog = 193:196; 
% nL2inj = unique(d2.cells(cL2inj,:));
% nL2erog = unique(d2.cells(cL2erog,:));
% cL4inj = 9409:9412;
% cL4erog = 193:196; 
% nL4inj = unique(d4.cells(cL4inj,:));
% % nL4erog = unique(d4.cells(cL4erog,:));
% 
% writeBCfiles('volErog_L2','VolumeForce','Dir','Flow',[],'erog_L2',[0 1 5 6 100],[0 -100 -100 0 0],4803);
% writeBCfiles('volErog_L4','VolumeForce','Dir','Flow',[],'erog_L4',[0 1 5 6 100],[0 -100 -100 0 0],4802);
% writeBCfiles('inj_L2_nod','NodeBC','Dir','Flow',[],'injection_L2',[0 1 5 6 100],[0 100 100 0 0],nL2inj);
% writeBCfiles('erog_L2_nod','NodeBC','Dir','Flow',[],'erog_L2',[0 1 5 6 100],[0 -100 -100 0 0],nL2erog);
% writeBCfiles('inj_L4_nod','NodeBC','Dir','Flow',[],'injection_L4',[0 1 5 6 100],[0 100 100 0 0],nL4inj);
% writeBCfiles('erog_L4_nod','NodeBC','Dir','Flow',[],'erog_L4',[0 1 5 6 100],[0 -100 -100 0 0],nL4erog);
% writeBCfiles('inj_L2_elem','ElementBC','Dir','Flow',[],'injection_L2',[0 1 5 6 100],[0 100 100 0 0],cL2inj);
% writeBCfiles('erog_L2_elem','ElementBC','Dir','Flow',[],'erog_L2',[0 1 5 6 100],[0 -100 -100 0 0],cL2erog);
% writeBCfiles('inj_L4_elem','ElementBC','Dir','Flow',[],'injection_L4',[0 1 5 6 100],[0 100 100 0 0],cL4inj);
% writeBCfiles('erog_L4_elem','ElementBC','Dir','Flow',[],'erog_L4',[0 1 5 6 100],[0 -100 -100 0 0],cL4erog);
% writeBCfiles('press_L4_nod','NodeBC','Dir','Flow',[],'press_L4',[0 1],[0 -100],[397,403]);
% writeBCfiles('press_L2_elem','ElementBC','Dir','Flow',[],'press_L2',[0 1],[0 -100],[9412,193]);
% writeBCfiles('press_L4_elem','ElementBC','Dir','Flow',[],'press_L4',[0 1],[0 -100],[9412,196]);
%writeBCfiles('bot_fix_L5','NodeBC','Dir','Poro',["x","y","z"],'fixed_bottom',0,0,bot_nodes);
%% Setup multidomain simulation
simParam = SimulationParameters('simParam.dat');
model = buildModelStruct('domains_mixed_FEMFV.dat',simParam);

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
model(1).Grid.topology.addSurface(1,'model-mixed/vol_1_surf_2.topol');
model(2).Grid.topology.addSurface(1,'model-mixed/vol_2_surf_1.topol');
model(2).Grid.topology.addSurface(2,'model-mixed/vol_2_surf_2.topol');
model(3).Grid.topology.addSurface(1,'model-mixed/vol_3_surf_1.topol');
model(3).Grid.topology.addSurface(2,'model-mixed/vol_3_surf_2.topol');
model(4).Grid.topology.addSurface(1,'model-mixed/vol_4_surf_1.topol');
model(4).Grid.topology.addSurface(2,'model-mixed/vol_4_surf_2.topol');
model(5).Grid.topology.addSurface(1,'model-mixed/vol_5_surf_1.topol');


%save('model.mat',"model")
%% Load mesh glue structure (avoid recomputing mortar every time)
fprintf('Computing mortar utils \n')
mG = MeshGlue(model,'interfaces_slaveFlow.dat');
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
%Solution loop
NSolv_MD.NonLinearLoop();

%%
for i = 1:numel(model) 
 mG.model(i).OutState.finalize()
end



