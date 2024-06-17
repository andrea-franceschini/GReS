clear
close all
%% contact search algorithm
% msh1 = Mesh();
% msh2 = Mesh();
% msh1.importGMSHmesh('3Dmesh/Mesh_coarseFlat.msh');
% msh2.importGMSHmesh('3Dmesh/Mesh_fineFlat.msh');


msh1 = Mesh();
msh2 = Mesh();
msh1.createCartesianGrid(2,[0 1],[0 1],10,10)
msh2.createCartesianGrid(2,[0 1],[0 1],14,14)
% get element connectivity between the teo meshes
cs = ContactSearching(msh1,msh2,18);

% analytical function on the master mesh
testFunc = @(x,y,z)  sin(3*x)+cos(3*y);
fIn = testFunc(msh1.coordinates(:,1), msh1.coordinates(:,2), msh1.coordinates(:,3));
fOutEx = testFunc(msh2.coordinates(:,1), msh2.coordinates(:,2), msh2.coordinates(:,3));
fID = fopen('out_nInt.dat', 'w');
nG = 15;
nInt = 15;
% compare with different numbers of GP
% compute Mortar operator
tIn = cputime;
[E,~,~] = compute_mortar3D(msh1, msh2, cs.elemConnectivity, nG, nInt);
t = cputime - tIn;
fOut = E*fIn;
errNorm = norm(fOut - fOutEx,2);


plotFunction(msh1, 'out_master', fIn)
plotFunction(msh2, 'out_slave', fOut)
plotFunction(msh2, 'out_slaveExact', fOutEx)


