%% pre processing and input data

clear
close all
clc

thisFolder = pwd;

% launch MRST
% cd('../../../MRST')
% addpath(genpath(fullfile(pwd,'testCornerPoint')))
%startup 
cd(thisFolder)

% grid parameters
gridName = 'cornerPointGrid';
dims = [28,44,28];
%refinement = [7,3];
nRockCells = 3;
%scale = [3e3,3e3,0.5e2];
scale = [1e2,1e2,10];

gresLog().setVerbosity(2);

%% process the grid
[mesh,newDims] = processCornerPointGrid(gridName,dims,nRockCells,scale);




%% run flow simulation

[mesh,pressures] = runFlowSimulation(mesh,newDims,nRockCells);

%% run mechanical simulation
%diary('log.txt')
%runMechanicsSimulation(mesh,pressures);
runContactMechanicsSimulation(mesh,pressures);
diary off


%% post processing

% displacement