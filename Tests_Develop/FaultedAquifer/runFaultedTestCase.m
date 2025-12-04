%% pre processing and input data
close all;
% clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% launch MRST
% cd('../../../MRST')
% addpath(genpath(fullfile(pwd,'testCornerPoint')))
%startup 


% grid parameters
gridName = 'Outputs/cornerPointGrid';
dims = [18,33,20];
nRockCells = 3;
scale = [1e2,1e2,10];

gresLog().setVerbosity(2);

%% process the grid
[mesh,newDims] = processCornerPointGrid(gridName,dims,nRockCells,scale);


%% run flow simulation
fileName = "flowCP.xml";
[mesh,pressures] = runFlowSimulation(fileName,mesh,newDims,nRockCells);


%% run mechanical simulation
fileName = "mechCP_stick.xml";
runContactMechanicsSimulation(fileName,mesh,pressures);
