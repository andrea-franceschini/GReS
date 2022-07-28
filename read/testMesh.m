close all;
clear;
clc;

% Setting the input file name
fileName = 'mesh.msh';

% Creation of an object of "Mesh"
mesh = Mesh();

% Starting a timer
tic;

% Calling a function from the class "Mesh" that import mesh data
mesh.importGMSHmesh(fileName);

% Reading the stopwatch timer
t1 = toc;

% Printing t1 = time to read 
fprintf('Time to read %.3f [s]\n', t1);

tic;
% Calling a function from the class "Mesh" for centroid calculation
mesh.finalize();
t2 = toc;

% Printing t2 = time to finalize
fprintf('Time to finalize %.3f [s]\n', t2);
list = fieldnames(mesh.surfaceRegions);
for i = 1 : length(list)
    fprintf('%s\n', cell2mat(list(i)));
end

