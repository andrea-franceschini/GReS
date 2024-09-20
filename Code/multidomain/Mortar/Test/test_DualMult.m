%%% TEST DUAL LAGRANGE MULTIPLIERS 

clear
close all
clc

type = 'gauss';
nGP = 3;
nInt = 4;

% interpolate variable from quad mesh to triangular mesh
nMaster = 3;
nSlave = 2;
x1 = 0;
x2 = 1;

master = zeros(nMaster,2);
slave = zeros(nSlave,2);

master(:,1) = (linspace(x1,x2,nMaster))';
slave(:,1) = (linspace(x1,x2,nSlave))';

master_top = build_topol(master(:,1));
slave_top = build_topol(slave(:,1));

mortar = Mortar2D(1,'set',master_top,slave_top,master,slave);
[D_RBF,M_RBF] = mortar.computeMortarRBF(nGP,nInt,type);
