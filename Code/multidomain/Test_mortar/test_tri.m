clear
close all
%%
% testing mortar class generalized to triangles and quadrilaterals

quadMaster = Mesh(); 
quadSlave = Mesh(); 
quadMaster.createCartesianGrid(2,1,[0 1],[0 1],8,8);
quadSlave.createCartesianGrid(2,1,[0 1],[0 1],5,5);
mQuad = Mortar3D(1,quadMaster,quadSlave);
[Dquad,Mquad,~,~,Equad] = mQuad.computeMortarRBF(2,4,'gauss');
f = @(x,y) sin(x) + cos(y);
fM = f(quadMaster.coordinates(:,1),quadMaster.coordinates(:,2));
fS = Equad*fM;
plotFunction(quadMaster,'outMQuad',fM)
plotFunction(quadSlave,'outSQuad',fS)

%%
f = @(x,y) sin(x) + cos(y);
triMaster = Mesh(); 
triSlave = Mesh();
triMaster.importGMSHmesh('triMaster.msh');
triSlave.importGMSHmesh('triSlave.msh');
mTri = Mortar3D(1,triMaster,1,triSlave,1);
[Dtri,Mtri,~,~,Etri] = mTri.computeMortarRBF(2,4,'gauss');
fM = f(triMaster.coordinates(:,1),triMaster.coordinates(:,2));
fS = Etri*fM;
plotFunction(triMaster,'outMTri',fM)
plotFunction(triSlave,'outSTri',fS)