%% pre processing and input data

clear
close all
clc

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);


%launch MRST
% cd('../../../MRST/mrst-2023a')
% addpath(genpath(fullfile(pwd,'testCornerPoint')))
% startup 
% cd(scriptDir)

% grid parameters
gridName = 'cornerPointGrid';
dims = [28,44,28];
%refinement = [7,3];
nRockCells = 3;
%scale = [3e3,3e3,0.5e2];
scale = [1e2,1e2,10];

%% process the grid
[mesh,newDims] = processCornerPointGrid(gridName,dims,nRockCells,scale);

% plot rescaled mesh for visualization purpose
% scaleNew = [5e2,5e2,1e2];
% meshVisual = mesh;
% meshVisual.coordinates = meshVisual.coordinates./scale;
% meshVisual.coordinates = meshVisual.coordinates.*scaleNew;




%% run flow simulation

[mesh,pressures] = runFlowSimulation(mesh,newDims,nRockCells);

%% run mechanical simulation
diary('log.txt')
%runMechanicsSimulation(mesh,pressures);
runContactMechanicsSimulation(mesh,pressures);
diary off


%% post processing

% displacement