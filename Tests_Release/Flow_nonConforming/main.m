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
meshFile = fullfile('Mesh','topLayer.vtk');
mesh = Mesh();
mesh.importMesh(meshFile);

% set up bc file for top layer
% get dirichlet node index (X = MIN, Y = MAX, Z = MAX)
tol = 1e-3;

[cx,cy,cz] = deal(mesh.coordinates(:,1),...
                  mesh.coordinates(:,2),...
                  mesh.coordinates(:,3));
                                          
X = min(cx);
Y = max(cy);
Z = max(cz);

dirNodTop= find(all(abs([cx-X, cy-Y cz-Z]) < tol,2));
writeBCfiles('BCs/bcTop','NodeBC','Dir','SinglePhaseFlow','topCorner_pressure',0,0,dirNodTop);


% read bottom layer mesh to set up boundary conditions
meshFile = fullfile('Mesh','bottomLayer.vtk');
mesh = Mesh();
mesh.importMesh(meshFile);

% set up bc file for bottom layer
% get dirichlet node index (X = MAX Y = MIN, Z = MIN)
[cx,cy,cz] = deal(mesh.coordinates(:,1),...
                  mesh.coordinates(:,2),...
                  mesh.coordinates(:,3));
                                          
X = max(cx);
Y = min(cy);
Z = min(cz);

dirNodBot= find(all(abs([cx-X, cy-Y cz-Z]) < tol,2));
writeBCfiles('BCs/bcBot','NodeBC','Dir','SinglePhaseFlow','topCorner_pressure',0,100,dirNodBot);

clear mesh


%% BUILD MODEL

% For now, simulation parameters are shared by all domains
simParam = SimulationParameters('simParam.dat');

% build model using domains input file (a shortcut to programmatically
% initialize separate model objects for each domain)
domainFile = fullfile('Domains','domain.xml');
domains = buildModelStruct(domainFile,simParam);

% Initialize the mortar utilities
interfFile = fullfile('Domains','interfaces.xml');
[interfaces,domains] = Mortar.buildInterfaceStruct(interfFile,domains);


%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(simParam,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();
