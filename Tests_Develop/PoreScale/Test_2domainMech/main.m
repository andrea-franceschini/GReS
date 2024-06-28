clear
close all

%% Testing multidomain model and inherent novel methods

fileName = "domains.dat";
% All the preprocessing is performed on different domains
% All istances of the classes are collected in a structure array
mod= buildModelStruct(fileName);

% SOLUTION ALGORITHM
fileName = 'testIntFile';
mG = MeshGlue(mod,fileName);

% get promechanics matrices from both domains
dt = 1; 
K1 = mod(1).Discretizer.getMat('Poro',mod(1).State);



