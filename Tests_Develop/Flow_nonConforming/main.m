clear
close all
%clc

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

%% PREPROCESSING TO WRITE BCS FILE PROGRAMMATICALLY

% read top layer mesh to set up boundary conditions
meshFile = fullfile('Mesh','Top_part.msh');
mesh = Mesh();
mesh.importMesh(meshFile);

elems = Elements(mesh,2);
% set up bc file for top layer
% get dirichlet surf index (X = 0, X = 590)
Xmin = 1e-3;
Xmax = 589.9;

[cx,cy,cz] = deal(elems.mesh.surfaceCentroid(:,1),...
                   elems.mesh.surfaceCentroid(:,2),...
                   elems.mesh.surfaceCentroid(:,3));

dirSurfX0 = find(cx<Xmin);
writeBCfiles('BCs/bcTop_H','SurfBC','Dir','SinglePhaseFlow','top_higher_pressure',0,50000,dirSurfX0);

dirSurfXM = find(cx>Xmax);
writeBCfiles('BCs/bcTop_L','SurfBC','Dir','SinglePhaseFlow','top_lower_pressure',0,45000,dirSurfXM);

% read bottom layer mesh to set up boundary conditions
meshFile = fullfile('Mesh','Bottom_part.msh');
mesh = Mesh();
mesh.importMesh(meshFile);

elems = Elements(mesh,2);
% set up bc file for bottom layer
% get dirichlet surf index (X = 0, X = 590)
Xmin = 1e-3;
Xmax = 589.9;

[cx,cy,cz] = deal(elems.mesh.surfaceCentroid(:,1),...
                  elems.mesh.surfaceCentroid(:,2),...
                  elems.mesh.surfaceCentroid(:,3));

dirSurfX0 = find(cx<Xmin);
writeBCfiles('BCs/bcBot_H','SurfBC','Dir','SinglePhaseFlow','Bottom_higher_pressure',0,50000,dirSurfX0);

dirSurfXM = find(cx>Xmax);
writeBCfiles('BCs/bcBot_L','SurfBC','Dir','SinglePhaseFlow','Bottom_lower_pressure',0,45000,dirSurfXM);

clear mesh

%% BUILD MODEL

% build model using domains input file (a shortcut to programmatically
% initialize separate model objects for each domain)
% domainFile = fullfile('Domains','domain.xml');
domainFile = fullfile('Domains','domain_new.xml');
domains = buildModel(domainFile);

domains(1).state.data.pressure(:) = 45000.;
domains(2).state.data.pressure(:) = 45000.;

% Initialize the mortar utilities
 % interfFile = fullfile('Domains','interfaces.xml');
interfFile = fullfile('Domains','interfaces_new.xml');
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);


%% RUN MODEL  
% A different solver is needed for models with non conforming domains
solver = MultidomainFCSolver(domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();





