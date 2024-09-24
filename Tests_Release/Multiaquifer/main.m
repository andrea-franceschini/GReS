% MULTI-AQUIFER MODEL: Non conforming grid blocks
%%
clear
close all

%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('MULTIAQUIFER MODEL \n')
fprintf('___________________\n\n')

%% Setup multidomain simulation
simParam = SimulationParameters('simParam.dat');
model = buildModelStruct('domains_mixed_FEMFV.dat',simParam);

model(1).Grid.topology.addSurface(1,'model-mixed/vol_1_surf_2.topol');
model(2).Grid.topology.addSurface(1,'model-mixed/vol_2_surf_1.topol');
model(2).Grid.topology.addSurface(2,'model-mixed/vol_2_surf_2.topol');
model(3).Grid.topology.addSurface(1,'model-mixed/vol_3_surf_1.topol');
model(3).Grid.topology.addSurface(2,'model-mixed/vol_3_surf_2.topol');
model(4).Grid.topology.addSurface(1,'model-mixed/vol_4_surf_1.topol');
model(4).Grid.topology.addSurface(2,'model-mixed/vol_4_surf_2.topol');
model(5).Grid.topology.addSurface(1,'model-mixed/vol_5_surf_1.topol');

%% Mortar interface computations
mG = MeshGlue(model,'interfaces_slaveFlow.dat');
fprintf('Done computing mortar utils\n')



%% modify interpolation operator to be sparse
for i = 1:numel(mG.interfaces)
  E = mG.interfaces(i).InterpOperator;
  E(abs(E) < 1e-4) = 0;
  E = E./sum(E,2);
  mG.interfaces(i).InterpOperator = sparse(E);
end

fprintf('Model Setup completed \n')
%% MULTIDOMAIN SOLUTION STRATEGY
NSolv_MD = NonLinearSolverMultiDomain(simParam,model,mG);
%Solution loop
NSolv_MD.NonLinearLoop();

%%
for i = 1:numel(model) 
 mG.model(i).OutState.finalize()
end



