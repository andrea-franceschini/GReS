%% pre processing and input data
clc
close all;
%clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir)
% launch MRST
% cd('../../../MRST')
% addpath(genpath(fullfile(pwd,'testCornerPoint')))
%startup 


% grid parameters
gridName = 'Outputs/cornerPointGrid';
dims = [24,40,26];
nRockCells = 3;
scale = [1e3,500,100];

gresLog().setVerbosity(3);

%% process the grid
[mesh,newDims] = processCornerPointGrid(gridName,dims,nRockCells,scale);


%% run flow simulation
fileName = "flowCP.xml";
[mesh,pressures] = runFlowSimulation(fileName,mesh,newDims,nRockCells);


%% run mechanical simulation
fileName = "mechCP.xml";
runContactMechanicsSimulation(fileName,mesh,pressures);
