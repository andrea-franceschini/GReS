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

% set up bc file for top layer
% get dirichlet node index (X = MIN, Y = MAX, Z = MAX)
% tol = 1e-3;
% 
% [cx,cy,cz] = deal(mesh.coordinates(:,1),...
%                   mesh.coordinates(:,2),...
%                   mesh.coordinates(:,3));
% 
% X = min(cx);
% Y = max(cy);
% Z = max(cz);
% 
% dirNodTop= find(all(abs([cx-X, cy-Y cz-Z]) < tol,2));
% % writeBCfiles('BCs/bcTop','NodeBC','Dir','SinglePhaseFlow','topCorner_pressure',0,0,dirNodTop);
% writeBCfiles('BCs/bcSup','NodeBC','Dir','SinglePhaseFlow','topCorner_pressure',0,0,dirNodTop);


% read bottom layer mesh to set up boundary conditions
%meshFile = fullfile('Mesh','bottomLayer.vtk');
meshFile = fullfile('Mesh','prova_inf.msh');
mesh = Mesh();
mesh.importMesh(meshFile);

% % set up bc file for bottom layer
% % get dirichlet node index (X = MAX Y = MIN, Z = MIN)
% [cx,cy,cz] = deal(mesh.coordinates(:,1),...
%                   mesh.coordinates(:,2),...
%                   mesh.coordinates(:,3));
% 
% X = max(cx);
% Y = min(cy);
% Z = min(cz);
% 
% dirNodBot= find(all(abs([cx-X, cy-Y cz-Z]) < tol,2));
% % writeBCfiles('BCs/bcBot','NodeBC','Dir','SinglePhaseFlow','topCorner_pressure',0,100,dirNodBot);
% writeBCfiles('BCs/bcInf','NodeBC','Dir','SinglePhaseFlow','topCorner_pressure',0,100,dirNodBot);

clear mesh


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
