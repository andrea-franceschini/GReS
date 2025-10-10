%% PORE SCALE SIMULATION
% Single phase flow with permeability depending on the channel size
% Pressure acting as a surface load on deformable grains
% Displacement of grains affecting the permeability change (loop)
% Mortar operators allow to transfer information from pressure to grains

% t = Mesh();
% t.importVTKmesh('mesh.vtk');
%%
clear
close all

% 
% t = Mesh();
% t.importMesh('Mesh/mesh.vtk');
%%
%[a,b,c] = mxImportGMSHmesh('Mesh/Fluid.msh');
f = Mesh();
f.importGMSHmesh('Mesh/Fluid_m.msh');

writeBCfiles('Fluid/dir_bot_flow','SurfBC','Dir','Flow',[],'botFix',0,5e3,f,1);
writeBCfiles('Fluid/dir_top_flow','SurfBC','Dir','Flow',[],'topFix',0,5.1e3,f,2);
% 
% 
g = Mesh();
g.importGMSHmesh('Mesh/Grain_m.msh'); 
% 
% % 
writeBCfiles('Grain/dir_fix_poroX','SurfBC','Dir','Poro',["x"],'fixX',0,0,g,1);
writeBCfiles('Grain/dir_fix_poroY','SurfBC','Dir','Poro',["y"],'fixY',0,0,g,2);
writeBCfiles('Grain/dir_fix_poroZ','SurfBC','Dir','Poro',["z"],'fixZ',0,0,g,3);

% set up simulation domains using multidomain class
model = buildModelStruct('domains.dat');

% set up interfaces for mortar computations from fluid to structure
glueF2S = MeshGlue(model,'interfaces_f2s.dat');
% set up interfaces for mortar computations from structure to fluid
glueS2F = MeshGlue(model,'interfaces_s2f.dat');


%% compute permeabilities for the flow domain
d = zeros(model(1).Grid.topology.nCells,1);
K = zeros(model(1).Grid.topology.nCells,1);
for i = 1:model(1).Grid.topology.nCells
   % compute channel size
   centroid = model(1).Grid.cells.cellCentroid(i,:);
   d(i) = computeChannelSize(centroid,glueF2S.interfaces(1).mortar.intMaster,glueF2S.interfaces(2).mortar.intMaster);
   K(i) = ((d(i))^3/12);
end
alpha = 0; 
poro = 0;

%% Solution algorithm
[rhsP,rhsU] = solvePoro(model,glueF2S.interfaces,glueS2F.interfaces,K,poro,alpha);

% Finilize print utilities
model(1).OutState.finalize();
model(2).OutState.finalize();

%% plot norm of displacements
u = reshape(model(2).State.dispCurr,3,[]);
u_norm = 1e6*(sqrt(sum(u.^2,1)));

plotFunction(g,'total_disp',u_norm')

%% convergence profiles on solutions at the first time step (ABSOLUTE RESIDUALS)


figure(1)
semilogy(1:length(rhsU),rhsU,'r',LineWidth=1,Marker='s',MarkerSize=8, MarkerFaceColor='r')
hold on
semilogy(1:length(rhsP)-1,rhsP(1:end-1),'b',LineWidth=1,Marker='^',MarkerSize=8, MarkerFaceColor='b')
xticks([1 2 3 4 5])
xlabel('iteration')
ylabel('Residual L^2 norm')
grid on

%% convergence profiles on solutions at the first time step (RELATIVE RESIDUALS)
figure(2)
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
semilogy(1:length(rhsU)-1,rhsU(1:end-1)./rhsU(1),'r',LineWidth=1,Marker='s',MarkerSize=8, MarkerFaceColor='r')
hold on
semilogy(1:length(rhsP)-2,rhsP(1:end-2)./rhsP(1),'b',LineWidth=1,Marker='^',MarkerSize=8, MarkerFaceColor='b')
xticks([1 2 3 4 5])
xlabel('Iteration')
ylabel('Relative residual $$L^2$$ norm')
l = legend('Mechanics','Flow');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'cmr10', 'FontSize', 12)
exportgraphics(gcf,'conv_profile2.pdf','ContentType','vector')