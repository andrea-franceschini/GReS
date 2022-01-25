close all;
clear;
clc;

fileName = 'mesh.msh';
mesh = Mesh();
tic;
mesh.importGMSHmesh(fileName);
t1 = toc;
fprintf('Time to read %.3f [s]\n', t1);
tic;
mesh.finalize();
t2 = toc;
fprintf('Time to finalize %.3f [s]\n', t2);
list = fieldnames(mesh.surfaceRegions);
for i = 1 : length(list)
    fprintf('%s\n', cell2mat(list(i)));
end
