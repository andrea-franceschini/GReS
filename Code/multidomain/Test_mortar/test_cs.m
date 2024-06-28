clear
close all
msh1 = Mesh();
msh2 = Mesh();

nM = 140;
nS = 200;

msh1.createCartesianGrid(2,1,[0 1],[0 1],nM,nM);
msh2.createCartesianGrid(2,1,[0 1],[0 1],nS,nS);
% Define object of 3D Mortar class
tic
mortar = Mortar3D(1,msh1,msh2);
t = toc;
plotFunction(msh1, 'test_cs_master', ones(msh1.nNodes,1));
plotFunction(msh2, 'test_cs_slave', ones(msh2.nNodes,1));