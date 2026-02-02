% close all;
% clear;
% 
% % Setting the input file name
% filevtk = 'meshVTK.vtk';
% 
% % Creation of an object of "Mesh"
% mesh = Mesh();
% 
% % Starting a timer
% tic;
% % Calling a function from the class "Mesh" that import mesh data
% %mesh.importGMSHmesh(filemsh);
% mesh.importVTKmesh(filevtk);
% t1 = toc;
% 
% % Printing t1 = time to read 
% fprintf('Time to read %.3f [s]\n', t1);
% 
% 
% % fprintf('Time to finalize %.3f [s]\n', t2);
% % list = fieldnames(mesh.surfaceRegions);
% % for i = 1 : length(list)
% %     fprintf('%s\n', cell2mat(list(i)));
% % end


clear 

m = structuredMesh(10,5,3,[0 50],[-10 40],[0 6]);


plotFunction(m,'testS1',ones(m.nNodes,1));

clear


xv = [0 5 15 20 25 27.5 30 31 32 38 43 50];
yv = [-10 0 10 20 25 30 40];
zv = [0 2 3 4 6];

m = structuredMesh(xv,yv,zv);

plotFunction(m,'testS2',ones(m.nNodes,1));


clear

b = BlockStructuredMesh(2,2,2,[0 1],[0 1],[0 1],3);
blockStructuedMesh = b.processGeometry();

plotFunction(blockStructuedMesh,'testB1',ones(blockStructuedMesh.nNodes,1));

clear


b = BlockStructuredMesh(2,2,2,[0 1],[0 1],[0 1],3);
b.refineRecursive([1 1 1],2);
b.refineRecursive([2 2 2],3);
b.refineRecursive([2 1 2],1);
blockStructuedMesh = b.processGeometry();

plotFunction(blockStructuedMesh,'testB2',ones(blockStructuedMesh.nNodes,1));



clear

b = BlockStructuredMesh(2,2,2,[0 1],[0 1],[0 1],3);
b.refineRecursive([1 1 1],2);
b.refineRecursive([2 2 2],3);
b.refineRecursive([2 1 2],1);
b.refineRecursive([2 1 2],2,2)
blockStructuedMesh = b.processGeometry();

plotFunction(blockStructuedMesh,'testB3',ones(blockStructuedMesh.nNodes,1));
