clear
close all
clc


scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);
%str_mod_ver = '2km'; % '2km' or '3km'


% write mesh files
elem_type = "hexa";
outFileTop = "topBlock";
outFileBottom = "bottomBlock";


[topMesh,bottomMesh] = deal(Mesh(),Mesh());

topMesh.importMesh(outFileTop + ".vtk");
bottomMesh.importMesh(outFileBottom + ".vtk");
% 
plotFunction(topMesh,'out_top',zeros(topMesh.nNodes,1));
plotFunction(bottomMesh,'out_bottom',zeros(bottomMesh.nNodes,1));

% write BC files
setBCfiles(topMesh,bottomMesh);

domainFile = 'Domains/domains.xml';
interfFile = 'Domains/interface.xml';
domains = buildModel(domainFile); 
% set verbosity 
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);


interfaces{1}.solvers(2).simparams.setVerbosity(2);

solver = ActiveSetContactSolver(domains,interfaces,5);

%solver.simParameters.setBackstepSkipFlag(1);

solver.NonLinearLoop();
solver.finalizeOutput();

%% plot profiles of multipliers along vertical axis (avoid opening paraview)
%plotStep(solver.results,2);




% 

