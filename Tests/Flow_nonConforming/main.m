% clear
close all
clc
input_dir = 'Inputs';

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
% cd(scriptDir);


%% BUILD MODEL

simparams = SimulationParameters(fullfile(input_dir,'flowNonConforming.xml'));
% Initialize the mortar utilities
[domains,interfaces] = buildModel(fullfile(input_dir,'flowNonConforming.xml'));


%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(simparams,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();
