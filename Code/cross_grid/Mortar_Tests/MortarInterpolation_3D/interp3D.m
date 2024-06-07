clear
close all

%type = 'gauss';
% msh1 = Mesh();
% msh2 = Mesh();
% msh1.importGMSHmesh('3Dmesh/Mesh_fineUnstruct.msh');
% msh2.importGMSHmesh('3Dmesh/Mesh_fineFlat.msh');
gauss = Gauss(12,3,3);
nM0 = 6;
nSizes = 5;
ratio = 0.6;
L2_eb = zeros(nSizes,1);
L2_rbf_w = L2_eb;
L2_rbf_g = L2_eb;

for i = 0:nSizes-1

    msh1 = Mesh();
    msh2 = Mesh();

    nM = nM0*2^i;
    nS = round(ratio*nM);

    msh1.createCartesianGrid(2,[0 1],[0 1],nM,nM);
    msh2.createCartesianGrid(2,[0 1],[0 1],nS,nS);
    % Define object of 3D Mortar class
    mortar = Mortar3D(1,msh1,msh2);
    %
    nG = 4;
    nInt = 4;

    [E_RBF_g,Mg,~,t1] = mortar.computeMortarRBF(nG,nInt,'gauss');
    [E_RBF_w,Mw,M1,t2] = mortar.computeMortarRBF(nG,nInt,'wendland');
    %[E2,M_eb,M2,tEB] = mortar.computeMortarElementBased(nG);

    % analytical function on the master mesh
    testFunc = @(x,y,z)  sin(3*x)+cos(3*y);
    fIn = testFunc(msh1.coordinates(:,1), msh1.coordinates(:,2), msh1.coordinates(:,3));
    fOutEx = testFunc(msh2.coordinates(:,1), msh2.coordinates(:,2), msh2.coordinates(:,3));
    fID = fopen('out_nInt.dat', 'w');
    fOutRBF_w = E_RBF_w*fIn;
    fOutRBF_g = E_RBF_g*fIn;
    %fOutEB = E2*fIn;
    % compute L2 interpolation error
    %L2_eb(i+1) = computeL2error(postProc(msh2,fOutEx,fOutEB,gauss));
    L2_rbf_w(i+1) = computeL2error(postProc(msh2,fOutEx,fOutRBF_w,gauss));
    L2_rbf_g(i+1) = computeL2error(postProc(msh2,fOutEx,fOutRBF_g,gauss));

end

%% plot convergence profiles
h = 1./(2.^(0:nSizes-1));
loglog(h,L2_rbf_g,'r-^')
hold on
loglog(h,L2_rbf_w,'b-^')
%loglog(h,L2_eb,'g-^')
legend('RBF - Rescaled Gauss', 'RBF - Rescaled Wendland','EB')
xlabel('mesh size')
ylabel('Quadratic error of interpolation')

plotFunction(msh1, 'out_master', fIn)
plotFunction(msh2, 'out_slaveRBF', fOutRBF_g)
plotFunction(msh2, 'out_slaveEB', fOutEB)



