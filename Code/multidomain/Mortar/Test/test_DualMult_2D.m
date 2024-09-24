%%% TEST DUAL LAGRANGE MULTIPLIERS for 3D problems
% provisional method computeMortarRBF_new
clear
close all

type = 'gauss';
gauss = Gauss(12,3,2);
nM0 = 2;
nSizes = 4;
ratio = 1.5;
L2_rbf_stand = zeros(nSizes,1);
L2_rbf_dual = L2_rbf_stand;
%L2_rbf_g = L2_eb;

for i = 0:nSizes-1

    msh1 = Mesh();
    msh2 = Mesh();

    nM = nM0*2^i;
    nS = round(ratio*nM);

    msh1.createCartesianGrid(2,1,[0 1],[0 1],nM,nM);
    msh2.createCartesianGrid(2,1,[0 1],[0 1],nS,nS);
    % Define object of 3D Mortar class
    mortar = Mortar3D(1,msh1,msh2);
    %
    nG = 3;
    nInt = 4;
    
    % analytical function on the master mesh
    testFunc = @(x,y,z)  sin(4*x).*cos(4*y);
    [Drbf_stand,Mrbf_stand] = mortar.computeMortarRBF_new(nG,nInt,type,'standard');
    [Drbf_dual,Mrbf_dual] = mortar.computeMortarRBF_new(nG,nInt,type,'dual');
    %[Deb,Meb] = mortar.computeMortarElementBased(nG);
    Erbf_stand =  Drbf_stand\Mrbf_stand;
    Erbf_dual =  Drbf_dual\Mrbf_dual;
    %fOutRBF = Erbf_stand*fIn;
    %fOutEB = Deb\(Meb*fIn);
    % compute L2 interpolation error
    %L2_eb(i+1) = computeL2error(postProc(msh2,fOutEx,fOutEB,gauss));
    %L2_rbf_w(i+1) = computeL2error(postProc(msh2,fOutEx,fOutRBF_w,gauss));
    L2_rbf_stand(i+1) = computeInterpError(mortar,Erbf_stand,testFunc);
    L2_rbf_dual(i+1) = computeInterpError(mortar,Erbf_dual,testFunc);
end

%% Save results in text file
% name = strcat('L2_',type,'_Int',num2str(nInt));
% fID = fopen(strcat('Results_lin\',name,'.dat'),'w');
% fprintf(fID,'%2.6e \n',L2_rbf);
% fID = fopen(strcat('Results_lin\L2_eb.dat'),'w');
% fprintf(fID,'%2.6e \n',L2_eb);
%% plot convergence profiles
h = 1./(ratio*nM0*2.^(0:nSizes-1));
fID = fopen('Results_lin\h.dat','w');
fprintf(fID,'%2.6e \n',h);
loglog(h,L2_rbf_stand,'r-^')
hold on
loglog(h,L2_dual,'g-^')
legend('standard','dual')
xlabel('mesh size')
ylabel('Quadratic error of interpolation')

%ratioW = L2_rbf_w(1:end-1)./L2_rbf_w(2:end);
%ratioEB = L2_eb(1:end-1)./L2_eb(2:end);
ratioRBF = L2_rbf(1:end-1)./L2_rbf(2:end);
% plotFunction(msh1, 'out_master', fIn)
% plotFunction(msh2, 'out_slaveRBF', fOutRBF_g)
% plotFunction(msh2, 'out_slaveEB', fOutEB)