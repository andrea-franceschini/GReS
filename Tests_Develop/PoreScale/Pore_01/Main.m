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

writeBCfiles('dir_bot_flow','SurfBC','Dir','Flow',[],'botFix',0,5e4,f,1);
writeBCfiles('dir_top_flow','SurfBC','Dir','Flow',[],'topFix',0,5.1e4,f,2);


g = Mesh();
g.importGMSHmesh('Mesh/Grain_m.msh'); 

% 
writeBCfiles('dir_fix_poro','SurfBC','Dir','Poro',["x"],'fix',0,0,g,1);
writeBCfiles('dir_fix_poroY','SurfBC','Dir','Poro',["y"],'fix',0,0,g,2);
writeBCfiles('dir_fix_poroZ','SurfBC','Dir','Poro',["z"],'fix',0,0,g,3);

% set up simulation domains using multidomain class
model = buildModelStruct('domains.dat');

% set up interfaces for mortar computations from fluid to structure
glueF2S = MeshGlue(model,'interfaces_f2s.dat');
% set up interfaces for mortar computations from strcture to fluid
glueS2F = MeshGlue(model,'interfaces_s2f.dat');


%% compute permeabilities for the flow domain
d = zeros(model(1).Grid.topology.nCells,1);
K = zeros(model(1).Grid.topology.nCells,1);
for i = 1:model(1).Grid.topology.nCells
   % compute channel size
   centroid = model(1).Grid.cells.cellCentroid(i,:);
   d(i) = computeChannelSize(centroid,glueF2S.interfaces(1).mortar.intMaster,glueF2S.interfaces(2).mortar.intMaster,30);
   K(i) = ((d(i))^3/12);
end
alpha = 0; 
poro = 0;

%% Solution algorithm
solvePoro(model,glueF2S.interfaces,glueS2F.interfaces,K,poro,alpha);

% Finilize print utilities
model(1).OutState.finalize();
model(2).OutState.finalize();


