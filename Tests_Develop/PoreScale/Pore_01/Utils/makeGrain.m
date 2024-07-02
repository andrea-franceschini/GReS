b = [0,0,0,100,100,100];
s = {[20,50,5,60];
   [90,50,110,60]};
% s = [20,50,5,60;
%    1,1,1,0.3;
%    0.85,0.4,0.5,0.3;
%    0,1,0.9,0.3];
genGrain('grainSphere.dat','fluidSphere.dat',b,s)

%%
% print meshes to paraview
g = Mesh();
f = Mesh();

g.importGMSHmesh('Grain_3.msh');
f.importGMSHmesh('Fluid_3.msh');
gfunc = ones(length(g.coordinates),1);
ffunc = ones(length(f.coordinates),1);

plotFunction(g,'outGrain',gfunc);
plotFunction(f,'outFluid',ffunc);