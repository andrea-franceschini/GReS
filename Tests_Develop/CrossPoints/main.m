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

% analytical solution
anal = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradx = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(2 - 2*X(1) + pi*cos(pi*X(1)));
grady = @(X) -pi*sin(pi*X(2)).*cos(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
gradz = @(X) -pi*cos(pi*X(2)).*sin(pi*X(3)).*(2*X(1)-X(1).^2 + sin(pi*X(1)));
h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
f = @(X) cos(pi*X(2)).*cos(pi*X(3)).*h(X(1));

% scaled solution along horizontal axis
% anal = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(4*X(1)-2*X(1).^2 + sin(2*pi*X(1)));
% gradx = @(X) cos(pi*X(2)).*cos(pi*X(3)).*(4 - 4*X(1) + 2*pi*cos(2*pi*X(1)));
% grady = @(X) -pi*sin(pi*X(2)).*cos(pi*X(3)).*(4*X(1)-2*X(1).^2 + sin(2*pi*X(1)));
% gradz = @(X) -pi*cos(pi*X(2)).*sin(pi*X(3)).*(4*X(1)-2*X(1).^2 + sin(2*pi*X(1)));
% h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
% f = @(X) cos(pi*X(2)).*cos(pi*X(3)).*h(X(1));

%% model commons
simParam = SimulationParameters('simParam.dat');

fname = 'cube';


meshFile = fullfile('Mesh','meshes',fname+".vtk");
mesh = Mesh();
mesh.importMesh(meshFile);

% set up bc file
nodes = unique(mesh.surfaces(ismember(mesh.surfaceTag,[9 10 11 12]),:));
c = mesh.coordinates(nodes,:);
vals = arrayfun(@(i) anal(c(i,:)),1:numel(nodes));
vals = reshape(vals,[],1);
writeBCfiles('BCs/bc','NodeBC','Dir','Poisson','manufactured_bc',0,0,nodes,vals);


% processing Poisson problem
domains = buildModelStruct('Domains/domain1block.xml',simParam);
domains.Discretizer.getSolver('Poisson').setAnalSolution(anal,f,gradx,grady,gradz);

[interfaces,domains] = Mortar.buildInterfaceStruct('Domains/interfaces_1.xml',domains);
% set up analytical solution

solver = MultidomainFCSolver(simParam,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();


