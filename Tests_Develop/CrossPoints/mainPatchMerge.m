clear
close all
clc

% this version of the cross-line patch test define the mortar surfaces with
% a L shape. This avoids creating copy of lagrange multipliers for the same
% node
% it serves to analyze the kernel of corner point non conformity

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);
% Change the current directory to the script's directory
cd(scriptDir);

%% model commons
simParam = SimulationParameters('simParam.dat');

fname = 'testPatchMerge2';


meshFile = fullfile('Mesh','meshes',fname+".vtk");
mesh = Mesh();
mesh.importMesh(meshFile);

% set up bc file
writeBCfiles('BCs/top_load','SurfBC','Neu',{'Poromechanics','z'},'top_load',0,-1,mesh,6);
writeBCfiles('BCs/bot_fix','NodeBC','Dir',{'Poromechanics','x','y','z'},'bot_fix',0,0,mesh,5);


% processing Poisson problem
domains = buildModelStruct('Domains/domain_patch_merge.xml',simParam);

[interfaces,domains] = Mortar.buildInterfaceStruct('Domains/interf_patch_merge.xml',domains);
% set up analytical solution

solver = MultidomainFCSolver(simParam,domains,interfaces);
solver.NonLinearLoop();
solver.finalizeOutput();


