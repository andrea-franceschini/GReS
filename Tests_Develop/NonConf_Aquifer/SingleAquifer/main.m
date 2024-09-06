% MULTI-AQUIFER MODEL: Non conforming grid blocks

%%
clear
close all

%% Define mesh objects
[top,res,bot] = deal(Mesh(),Mesh(),Mesh());
top.importGMSHmesh('Mesh/Overburden.msh')
res.importGMSHmesh('Mesh/Reservoir.msh')
bot.importGMSHmesh('Mesh/Underburden.msh')
%% Utils
%writeBCfiles('bottom_fixed','SurfBC','Dir','Poro',["x","y","z"],'botFix',0,0,bot,2);
% writeBCfiles('BCs/press_inj','NodeBC','Dir','Flow',[],'inj',0,10,421);
% writeBCfiles('BCs/press_erog','NodeBC','Dir','Flow',[],'erog',0,-10,574);
%writeBCfiles('point_erog','NodeBC','Dir','Flow',[],'dir_flow',[0 1 5 100],[0 -100 -500 -500],3328);
%writeBCfiles('flux','VolumeForce',[],'Flow',[],'flux_bc',[0 1 5 100],[0 -0.1 -0.1 -0.1],3388);
% p = find(all([res.coordinates(:,1)==1000,res.coordinates(:,2)==1000],2));
% p_c = res.coordinates(p,:);

%% Setup multidomain simulation
simParam = SimulationParameters('simParam.dat');
model = buildModelStruct('domains.dat',simParam);

%% Setup computations between non conforming domains
mG = MeshGlue(model,'interfaces.dat');

%% MULTIDOMAIN SOLUTION STRATEGY
NSolv_MD = NonLinearSolverMultiDomain(simParam,model,mG);
% Solution loop
NSolv_MD.NonLinearLoop();

for i = 1:numel(model)
   model(i).OutState.finalize()
end


