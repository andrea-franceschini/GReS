% clear
close all
clc
input_dir = 'Input';

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
% cd(scriptDir);


%% BUILD MODEL
fName = fullfile(input_dir,'flowNonConforming.xml');

simparams = SimulationParameters(fName);
printUtils = OutState(fName);
% Initialize the mortar utilities
[domains,interfaces] = buildModel(fName);


%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = NonLinearImplicit('simulationparameters',simparams,...
                           'domains',domains,...
                           'interface',interfaces, ...
                           'output',printUtils);
solver.simulationLoop();
