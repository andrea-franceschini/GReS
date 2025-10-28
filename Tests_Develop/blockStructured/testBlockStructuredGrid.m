clear
close all
clc

% block structured mesh object 
% INPUT: 
% the limit of the box grid
% the number of macro-block on each direction
% the maximum refinement level
% Refinement is achieved with an octree strategy
fprintf('Processing the block structured mesh')

% create the base conforming structured mesh
blockMesh = BlockStructuredMesh([0,1;0,1;0,0.25],[3,3,1],4);

% OPTION 1: refine recursively all children of target macro block
blockMesh.refineRecursive([2,2,1],2);
% blockMesh.refineRecursive([2,2,1],3,3);
% blockMesh.refineRecursive([2,2,1],3,7);
% blockMesh.refineRecursive([1,2,1],3);
% blockMesh.refineRecursive([3,1,1],2);
% blockMesh.refineRecursive([3,3,1],1);

% examples
% case 1
%blockMesh.refineRecursive([2,2,1],1);
% case 1
%blockMesh.refineRecursive([2,2,1],2);
% case 3
%blockMesh.refineRecursive([2,2,1],2);
% blockMesh.refineRecursive([2,2,1],3,3);
% blockMesh.refineRecursive([2,2,1],3,7);


% OPTION 2: refine around a target point until the maximum refinement level
% is reached
% blockMesh.refineRecursive([0.5,0.5,0.1],1); 
% blockMesh.refineRecursive([0.4,0.6,0.1],3); 
% blockMesh.refineGrid([0.3,0.6,0.06]);


% process a GReS mesh object based on the block structured grid
% Mortar surfaces are automatically generated
% surfacTag 1 -> master
% surfaceTag 2 -> slave
mesh = blockMesh.processGeometry();

% plot the mesh object for testing
plotFunction(mesh,'test',ones(mesh.nNodes,1));

meshMaster = getSurfaceMesh(mesh,1);
meshSlave = getSurfaceMesh(mesh,2);
plotFunction(meshSlave,'testSlave',ones(meshSlave.nNodes,1));
plotFunction(meshMaster,'testMaster',ones(meshMaster.nNodes,1));

% define model 
model = ModelType("Poromechanics_FEM");
simParam = SimulationParameters("simParam.dat");
elems = Elements(mesh,2);
faces = Faces(model, mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);
mat = Materials(model,"materialsList.dat");
dofmanager = DoFManager(mesh,model);

% Create and set the print utility
printUtils = OutState(model,mesh,'outTime.dat','folderName','Output_multiBlock','flagMatFile',false);


% Write BC files programmatically
mkdir BCs
foldName = 'BCs';
writeBCfiles(strcat(foldName,'/topLoad'),'SurfBC','Neu',{'Poromechanics','y'},'TopLoad',0,-1,mesh,6);
writeBCfiles(strcat(foldName,'/fixBot'),'NodeBC','Dir',{'Poromechanics','x','y','z'},'fixedBot',0,0,mesh,5);
writeBCfiles(strcat(foldName,'/fixZ'),'NodeBC','Dir',{'Poromechanics','z'},'fixedZ',0,0,mesh,[3,4]);
writeBCfiles(strcat(foldName,'/fixX'),'NodeBC','Dir',{'Poromechanics','x'},'fixedX',0,0,mesh,[7,8]);

% Collect BC input file in a list
fileName = ["BCs/topLoad.dat","BCs/fixBot.dat"];

%
% Create an object of the "Boundaries" class 
bound = Boundaries(fileName,model,grid);


% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

% build the mortar interface with xml input
interfFile = fullfile('Domains','interfaces.xml');
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domain);


% Print model initial state
printState(domain);

% Istanciate a mortar solver
Solver = MultidomainFCSolver(domain,interfaces);

% Solve the problem
Solver.NonLinearLoop();

% Finalize the print utility
Solver.finalizeOutput();




