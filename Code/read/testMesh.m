close all;
clear;

% Setting the input file name
filevtk = 'BottomBlock_tetra.vtk';

% Creation of an object of "Mesh"
mesh = Mesh();

% Starting a timer
tic;
% Calling a function from the class "Mesh" that import mesh data
%mesh.importGMSHmesh(filemsh);
mesh.importVTKmesh(filevtk);
t1 = toc;

% Printing t1 = time to read 
fprintf('Time to read %.3f [s]\n', t1);


% fprintf('Time to finalize %.3f [s]\n', t2);
% list = fieldnames(mesh.surfaceRegions);
% for i = 1 : length(list)
%     fprintf('%s\n', cell2mat(list(i)));
% end

