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
% writeBCfiles('BCs/bottom_fixed','SurfBC','Dir','Poro',["x","y","z"],'botFix',0,0,bot,1);
% writeBCfiles('BCs/press_inj','NodeBC','Dir','Flow',[],'inj',0,10,421);
% writeBCfiles('BCs/press_erog','NodeBC','Dir','Flow',[],'erog',0,-10,574);
% p = find(all([res.coordinates(:,1)==1000,res.coordinates(:,2)==1000],2));
% p_c = res.coordinates(p,:);

%% Setup multidomain simulation
simParam = SimulationParameters()
model = buildModelStruct('domains.dat');

%% Setup computations between non conforming domains
mG = MeshGlue(model,'interfaces.dat');

%% MULTIDOMAIN SOLUTION STRATEGY
NSolv_MD = NonLinearSolverMultiDomain(model,mG.interfaces);
% Solution loop
NSolv_MD.NonLinearLoop();


