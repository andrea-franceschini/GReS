clear
close all
clc

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

%% PREPROCESSING TO WRITE BCS FILE PROGRAMMATICALLY

% read top layer mesh to set up boundary conditions
% meshFile = fullfile('Mesh','topLayer.vtk');
 meshFile = fullfile('Mesh','prova_sup.msh');
mesh = Mesh();
mesh.importMesh(meshFile);

% read bottom layer mesh to set up boundary conditions
%meshFile = fullfile('Mesh','bottomLayer.vtk');
meshFile = fullfile('Mesh','prova_inf.msh');
mesh = Mesh();
mesh.importMesh(meshFile);


clear mesh

domains(1).state.data.pressure(:) = 0.;
domains(2).state.data.pressure(:) = 0.;

%% BUILD MODEL

% build model using domains input file (a shortcut to programmatically
% initialize separate model objects for each domain)
% domainFile = fullfile('Domains','domain.xml');
domainFile = fullfile('Domains','domain_new.xml');
domains = buildModel(domainFile);

% Initialize the mortar utilities
 % interfFile = fullfile('Domains','interfaces.xml');
interfFile = fullfile('Domains','interfaces_new.xml');
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);


%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();
