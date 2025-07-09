clear
close all
clc

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

%% Poisson problem with single domain in 3D. Testing new poisson module

% analytical solution of Poisson problem
anal = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradx = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2 - 2*X(1) + pi*cos(pi*X(1)));
grady = @(X) -pi*sin(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradz = @(X) -pi*cos(pi*X(2)).*sin(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
f = @(X) cos(pi*X(2)).*cos(pi*X(3)).*h(X(1));

simParam = SimulationParameters('simParam.dat');


meshFile = fullfile('Mesh','domain.vtk');
mesh = Mesh();
mesh.importMesh(meshFile);

% set up bc file
nodes = unique(mesh.surfaces(ismember(mesh.surfaceTag,[2 4]),:));
c = mesh.coordinates(nodes,:);
vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
vals = reshape(vals,[],1);
writeBCfiles('BCs/bc','NodeBC','Dir','Poisson','manufactured_bc',0,0,nodes,vals);

clear mesh

% build model and set up Poisson manufactured solution
domainFile = fullfile('Domains','domain.xml');
domains = buildModelStruct(domainFile,simParam);
domains.Discretizer.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);

% build interface structure
interfFile = fullfile('Domains','interfaces.xml');
[interfaces,domains] = Mortar.buildInterfaceStruct(interfFile,domains);

solver = MultidomainFCSolver(simParam,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();
