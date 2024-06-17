clear
close all

type = 'gauss';
msh1 = Mesh();
msh2 = Mesh();
msh1.importGMSHmesh('3Dmesh/Mesh_coarse_HEXA.msh');
msh2.importGMSHmesh('3Dmesh/Mesh_fine_HEXA.msh');
gauss = Gauss(12,3,3);
% nM0 = 8;
% nSizes = 6;
% ratio = 0.4;
% L2_eb = zeros(nSizes,1);
% L2_rbf_w = L2_eb;
% L2_rbf_g = L2_eb;

%Define object of 3D Mortar class
mortar = Mortar3D(1,msh1,1,msh2,1);

nG = 6;
nInt = 4;

[E_RBF_g,Mg,~,tg] = mortar.computeMortarRBF(nG,nInt,'gauss');
[E_RBF_w,Mw,~,tw] = mortar.computeMortarRBF(nG,nInt,'wendland');
[E2,Meb,M2,t2] = mortar.computeMortarElementBased(nG);

%analytical function on the master mesh
testFunc = @(x,y,z)  sin(3*x)+cos(3*x);
fIn = testFunc(msh1.coordinates(:,1), msh1.coordinates(:,2), msh1.coordinates(:,3));
fOutEx = testFunc(msh2.coordinates(:,1), msh2.coordinates(:,2), msh2.coordinates(:,3));
% fID = fopen('out_nInt.dat', 'w');
fOutRBF_w = E_RBF_w*fIn;
fOutRBF_g = E_RBF_g*fIn;
fOutEB = E2*fIn;
% compute L2 interpolation error
L2_eb = computeL2error(postProc(msh2,fOutEx,fOutEB,gauss));
L2_rbf_w = computeL2error(postProc(msh2,fOutEx,fOutRBF_w,gauss));
L2_rbf_g = computeL2error(postProc(msh2,fOutEx,fOutRBF_g,gauss));


plotFunction(msh1, 'out_master', fIn)
plotFunction(msh2, 'out_slaveRBF_w', fOutRBF_w)
plotFunction(msh2, 'out_slaveRBF_g', fOutRBF_g)
plotFunction(msh2, 'out_slaveEB', fOutEB)
plotFunction(msh2, 'out_slaveExact', fOutEx)
plotFunction(msh2, 'out_RelErr', rel_g)


