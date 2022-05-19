close all;
clear;
clc;

% Setting input file
fileName = 'boundaries.dat';

%Creation of an object of "Boundaries"
bound = Boundaries(fileName);
% Starting a timer
tic;

disp = bound.getBC('dir');
disp.DirBoundary()
% nodeforce = bound.getBC('neumann');
% nodeforce.NeuBoundary()

% Reading the stopwatch timer
t1 = toc;
% Printing t1 = time to read 
fprintf('Time to read %.3f [s]\n', t1);
