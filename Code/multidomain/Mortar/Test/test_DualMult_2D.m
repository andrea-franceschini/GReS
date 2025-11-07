%%% TEST DUAL LAGRANGE MULTIPLIERS for 3D problems
% provisional method computeMortarRBF_new
clear
close all
clc

type = 'gauss';
gauss = Gauss(12,3,2);
nM0 = 4;
nSizes = 3;
ratio = 0.6;
L2_rbf_stand = zeros(nSizes,1);
L2_rbf_dual = L2_rbf_stand;
L2_eb = L2_rbf_stand;
%L2_rbf_g = L2_eb;

for i = 0:nSizes-1

    msh1 = Mesh();
    msh2 = Mesh();

    nM = nM0*2^i;
    nS = round(ratio*nM);

    msh1.createCartesianGrid(2,1,[0 1],[0 1],nM,nM);
    msh2.createCartesianGrid(2,1,[0.1 1.01],[0.1 1.05],nS,nS);
    % Define object of 3D Mortar class
    mortar = Mortar3D(1,msh1,msh2);
    %
    nG = 4;
    nInt = 4;

    % analytical function on the master mesh
    testFunc = @(x,y,z)  sin(4*x).*cos(4*y);
    fprintf('Standard meshRef %i \n',i)
    [Drbf_stand,Mrbf_stand] = mortar.computeMortarRBF(nG,nInt,type,'standard');
    fprintf('________________________________ \n');
    fprintf('Dual meshRef %i \n',i)
    [Drbf_dual,Mrbf_dual] = mortar.computeMortarRBF(nG,nInt,type,'dual');
    fprintf('________________________________ \n');
    fprintf('EB meshRef %i \n',i)
    [D_eb,M_eb] = mortar.computeMortarElementBased(nG);
    fprintf('________________________________ \n');
    % Dtest = full(Drbf_dual);
    [Deb,Meb] = mortar.computeMortarElementBased(nG);
    Erbf_stand =  Drbf_stand\Mrbf_stand;
    Erbf_dual =  Drbf_dual\Mrbf_dual;
    E_eb = D_eb\M_eb;
    diag_meas = norm(Drbf_dual-diag(diag(Drbf_dual)),'fro');
    isPU_stand = norm(abs(sum(Erbf_stand,2)-ones(size(Erbf_stand,1),1)),2);
    isPU_dual = norm(abs(sum(Erbf_dual,2)-ones(size(Erbf_dual,1),1)),2);

    fIn = testFunc(msh1.coordinates(:,1),msh1.coordinates(:,2),msh1.coordinates(:,3));
    fOutEx = testFunc(msh2.coordinates(:,1),msh2.coordinates(:,2),msh2.coordinates(:,3));
    fOutstand = Erbf_stand*fIn;
    fOutdual = Erbf_dual*fIn;
    fOutEB = E_eb*fIn;
    L2_rbf_stand(i+1) = computeL2error(postProc(msh2,fOutEx,fOutstand,gauss));
    L2_rbf_dual(i+1) = computeL2error(postProc(msh2,fOutEx,fOutdual,gauss));
    L2_eb(i+1) = computeL2error(postProc(msh2,fOutEx,fOutEB,gauss));
end

%% Save results in text file
% name = strcat('L2_',type,'_Int',num2str(nInt));
% fID = fopen(strcat('Results_lin\',name,'.dat'),'w');
% fprintf(fID,'%2.6e \n',L2_rbf);
% fID = fopen(strcat('Results_lin\L2_eb.dat'),'w');
% fprintf(fID,'%2.6e \n',L2_eb);
%% plot convergence profiles
h = 1./(ratio*nM0*2.^(0:nSizes-1));
fID = fopen('Results_lin_h.dat','w');
fprintf(fID,'%2.6e \n',h);
loglog(h,L2_rbf_stand,'r-^')
hold on
loglog(h,L2_rbf_dual,'g-^')
loglog(h,L2_eb,'k-s')
legend('standard','dual','EB')
xlabel('mesh size')
ylabel('Quadratic error of interpolation')

%ratioW = L2_rbf_w(1:end-1)./L2_rbf_w(2:end);
%ratioEB = L2_eb(1:end-1)./L2_eb(2:end);
ratioRBF_stand = L2_rbf_stand(1:end-1)./L2_rbf_stand(2:end);
ratioRBF_dual = L2_rbf_dual(1:end-1)./L2_rbf_dual(2:end);
% plotFunction(msh1, 'out_master', fIn)
% plotFunction(msh2, 'out_slaveRBF', fOutRBF_g)
% plotFunction(msh2, 'out_slaveEB', fOutEB)