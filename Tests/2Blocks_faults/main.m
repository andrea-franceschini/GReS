clear
close all

% fault material parameters
coes = 0;
phi = 30; % degrees
 
scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);
%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('2blocks model \n')
fprintf('___________________\n\n')
[leftMesh,rightMesh] = deal(Mesh(),Mesh());


meshLeftFile = 'Mesh/LeftBlock_hexa.msh';
meshRightFile = 'Mesh/RightBlock_hexa.msh';
%domainFile = 'Domains/domains_hexa_P0.dat';

% run geo files
mshLeft = 'Mesh/LeftBlockHexa.geo';
mshRight = 'Mesh/RightBlockHexa.geo';
cmd = sprintf('gmsh - %s',mshLeft);
[status, result] = system(cmd);
cmd = sprintf('gmsh - %s',mshRight);
[status, result] = system(cmd);

% mshTest = Mesh();
% mshTest.importMesh('Mesh/RightBlockHexa.vtk');

leftMesh.importGMSHmesh(meshLeftFile);
rightMesh.importGMSHmesh(meshRightFile);
% 
% plotFunction(leftMesh,'out_L',zeros(leftMesh.nNodes,1),"sol");
% plotFunction(rightMesh,'out_R',zeros(rightMesh.nNodes,1),"sol");

% write BC files
setBCfiles(leftMesh,rightMesh);

domainFile = 'Domains/domains.xml';
interfFile = 'interface.xml';
simParam = SimulationParameters('simParam.dat');
domains = buildModelStruct(domainFile,simParam);
%scale domain size to millimiters
[interfaces,domains] = Mortar.buildInterfaceStruct(interfFile,domains);

% setting initial multiplier value manually
% compressive sigma_n = 1 kPa for entire interface
interfaces{1}.iniMultipliers{1}(1:3:end) = -1; 
interfaces{1}.multipliers(1).curr(1:3:end) = -1; 
interfaces{1}.multipliers(1).prev(1:3:end) = -1; 

solver = MultidomainFCSolver(simParam,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();

%% plot profiles of multipliers along vertical axis (avoid opening paraview)
plotStep(solver.results,2);




% 

