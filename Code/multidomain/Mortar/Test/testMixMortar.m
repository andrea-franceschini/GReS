clear
close all
clc

% interpolate variable from quad mesh to triangular mesh
tri = Mesh();
quad = Mesh();

quad.createCartesianGrid(2,1,[0 1],[0 1],7,7);
tri.importGMSHmesh('tri.msh');
m = Mortar3D(1,quad,1,tri,1);
[D,M,~,~,E] = m.computeMortarRBF_new(2,4,'gauss','dual');
f = @(x,y) sin(x)+cos(y);
fM = f(quad.coordinates(:,1),quad.coordinates(:,2));
fS = E*fM;
plotFunction(quad,'outQuad',fM);
plotFunction(tri,'outTri',fS);