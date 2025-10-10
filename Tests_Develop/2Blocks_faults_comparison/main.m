clear
close all
clc

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

% write mesh files
elem_type = "hexa";
outFileLeft = "leftBlock";
outFileRight = "rightBlock";
NXL = 3; NYL = 10; NZL = 12; 
NXR = 4; NYR = 10; NZR = 12; 

comLeft = "python Mesh/leftBlock.py " + outFileLeft + " " + num2str(NXL) + " " + num2str(NYL) + " "  + num2str(NZL) + " " + elem_type;
comRight = "python Mesh/rightBlock.py " + outFileRight + " " + num2str(NXR) + " " + num2str(NYR) + " "  + num2str(NZR) + " " + elem_type;

system(comLeft);
system(comRight);

leftMesh.importMesh(outFileLeft + ".vtk");
rightMesh.importMesh(outFileRight + ".vtk");
% 
% plotFunction(leftMesh,'out_L',zeros(leftMesh.nNodes,1),"sol");
% plotFunction(rightMesh,'out_R',zeros(rightMesh.nNodes,1),"sol");

% write BC files
setBCfiles(leftMesh,rightMesh);

domainFile = 'Domains/domains.xml';
interfFile = 'Domains/interface_stick.xml';
domains = buildModel(domainFile); 
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);

% setting initial normal traction
% compressive sigma_n = 1 kPa for entire interface
interfaces{1}.iniMultipliers(1:3:end) = -1; 
interfaces{1}.multipliers.curr(1:3:end) = -1; 
interfaces{1}.multipliers.prev(1:3:end) = -1; 

% set verbosity 
domains(2).simparams.setVerbosity(2);



solver = MultidomainFCSolver(domains,interfaces);


solver.NonLinearLoop();
solver.finalizeOutput();

%% plot profiles of multipliers along vertical axis (avoid opening paraview)
%plotStep(solver.results,2);


% building L2 projection operator to try smooth out the multipliers
% simple approach (weighed area)




% 

