%%% TEST DUAL LAGRANGE MULTIPLIERS 

clear
close all
clc

type = 'gauss';
nGP = 3;
nInt = 4;
f = @(x) sin(3*x);

% interpolate variable from quad mesh to triangular mesh
nMaster0 = 2;
r = 3;
ndim = 7;
err_standard = zeros(ndim,1);
err_dual = zeros(ndim,1);
hs = err_dual;

x1 = 0;
x2 = 1;
for i = 1:ndim
    nMaster = nMaster0*2^i;
    nSlave = round(r*nMaster);
    master = zeros(nMaster,2);
    slave = zeros(nSlave,2);
    hs(i) = (x2-x1)/nSlave;

    master(:,1) = (linspace(x1,x2,nMaster))';
    slave(:,1) = (linspace(x1,x2,nSlave))';

    master_top = build_topol(master(:,1));
    slave_top = build_topol(slave(:,1));

    mortar = Mortar2D(1,'set',master_top,slave_top,master,slave);
    tic
    [D_RBF_standard,M_RBF_standard] = mortar.computeMortarRBF(nGP,nInt,type,'standard');
    t_st = toc;
    E_RBF_standard = D_RBF_standard\M_RBF_standard;
    tic
    [D_RBF_dual,M_RBF_dual] = mortar.computeMortarRBF(nGP,nInt,type,'dual');
    t_dual = toc;
    E_RBF_dual = D_RBF_dual\M_RBF_dual;
    err_standard(i) = computeInterpError(mortar,E_RBF_standard,f);
    err_dual(i) = computeInterpError(mortar,E_RBF_dual,f);
end

loglog(hs,err_standard,'r-o')
hold on
loglog(hs,err_dual,'g--s')
legend('standard','dual')
rat_stand = err_standard(1:end-1)./err_standard(2:end);




