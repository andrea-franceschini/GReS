clear
close all


type = 'gauss';
msh1 = Mesh();
msh2 = Mesh();
nM = 6;
nS = 14;
msh1.createCartesianGrid(2,1,[-1 1],[-1 1],nM,nM);

%plotFunction(msh1, 'test', ones(msh1.nNodes,1))

msh2.createCartesianGrid(2,1,[-1 1],[-1 1],nS,nS);

% change the shape of the surfaces with z = f(x,y)
k=0.5;
shape = @(x,y) 1-k*x.^2-k*y.^2;
msh1.coordinates(:,3) = shape(msh1.coordinates(:,1),msh1.coordinates(:,2));
msh2.coordinates(:,3) = shape(msh2.coordinates(:,1),msh2.coordinates(:,2));

% msh2.coordinates(:,3) = 0.25;
gauss = Gauss(12,3,2);
% nM0 = 8;
% nSizes = 6;
% ratio = 0.4;
% L2_eb = zeros(nSizes,1);
% L2_rbf_w = L2_eb;
% L2_rbf_g = L2_eb;

%Define object of 3D Mortar class
mortar = Mortar3D(1,msh1,msh2);

nG = 2;
nInt = 4;

[D_rbf,M_rbf] = mortar.computeMortarRBF(nG,nInt,'gauss');
%[E_RBF_w,Mw,~,tw] = mortar.computeMortarRBF(nG,nInt,'wendland');
[D_eb,M_eb] = mortar.computeMortarElementBased(nG);

%analytical function on the master mesh
testFunc = @(x,y,z)  sin(x)+cos(y);
fIn = testFunc(msh1.coordinates(:,1), msh1.coordinates(:,2), msh1.coordinates(:,3));
fOutEx = testFunc(msh2.coordinates(:,1), msh2.coordinates(:,2), msh2.coordinates(:,3));
% fID = fopen('out_nInt.dat', 'w');
%fOutRBF_w = E_RBF_w*fIn;
fOutRBF_g = D_rbf\(M_rbf*fIn);
fOutEB = D_eb\(M_eb*fIn);
% compute L2 interpolation error
L2_eb = computeL2error(postProc(msh2,fOutEx,fOutEB,gauss));
%L2_rbf_w = computeL2error(postProc(msh2,fOutEx,fOutRBF_w,gauss));
L2_rbf_g = computeL2error(postProc(msh2,fOutEx,fOutRBF_g,gauss));

rel = 100*abs(L2_eb-L2_rbf_g)/L2_eb;
plotFunction(msh1, 'out_master', fIn)
% plotFunction(msh2, 'out_slaveRBF_w', fOutRBF_w)
plotFunction(msh2, 'out_slaveRBF_g', fOutRBF_g)
plotFunction(msh2, 'out_slaveEB', fOutEB)
%plotFunction(msh2, 'out_slaveRBF_g', fOutRBF_g)





