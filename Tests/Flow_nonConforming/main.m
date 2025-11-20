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
% build model using domains input file (a shortcut to programmatically
% initialize separate model objects for each domain)
domains = buildModel(fullfile(input_dir,'domain.xml'));

% Initialize the mortar utilities
[interfaces,domains] = Mortar.buildInterfaces(fullfile(input_dir,'interfaces.xml'),domains);


%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();
