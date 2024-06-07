clear
close all

type = 'gauss';
msh1 = Mesh();
msh2 = Mesh();
nM = 8;
nS = 12;
msh1.createCartesianGrid(2,[-1 1],[-1 1],nM,nM);
msh2.createCartesianGrid(2,[-1 1],[-1 1],nS,nS);

% change the shape of the surfaces with z = f(x,y)
% k=0.5;
% shape = @(x,y) 1-k*x.^2-k*y.^2;
% msh1.coordinates(:,3) = shape(msh1.coordinates(:,1),msh1.coordinates(:,2));
% msh2.coordinates(:,3) = 0.12+shape(msh2.coordinates(:,1),msh2.coordinates(:,2));

msh2.coordinates(:,3) = 0.25;
gauss = Gauss(12,3,3);
% nM0 = 8;
% nSizes = 6;
% ratio = 0.4;
% L2_eb = zeros(nSizes,1);
% L2_rbf_w = L2_eb;
% L2_rbf_g = L2_eb;

%Define object of 3D Mortar class
mortar = Mortar3D(1,msh1,msh2);

nG = 4;
nInt = 6;

[E_RBF_g,Mg,~,tg] = mortar.computeMortarRBF(nG,nInt,'gauss');
[E_RBF_w,Mw,~,tw] = mortar.computeMortarRBF(nG,nInt,'wendland');
[E2,Meb,M2,t2] = mortar.computeMortarElementBased(nG);

%analytical function on the master mesh
testFunc = @(x,y,z)  sin(3*x)+cos(3*y);
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
plotFunction(msh2, 'out_slave_eb', fOutEB)
plotFunction(msh2, 'out_slaveExact', fOutEx)
%plotFunction(msh2, 'out_RelErr', rel_g)



