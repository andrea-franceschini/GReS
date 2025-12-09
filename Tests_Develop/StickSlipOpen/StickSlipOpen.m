clc
close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir)


% set mesh 
b1 = BlockStructuredMesh([0,4;0 20;0 20],[4 10 10],1);
meshL = processGeometry(b1);

b2 = BlockStructuredMesh([4,8;0 20;0 20],[4 15 15],1);
meshR = processGeometry(b2);

% define model 

simParam = SimulationParameters("StickSlipOpen.xml");

elemsL = Elements(meshL,2);
facesL = Faces(meshL);
gridL = struct('topology',meshL,'cells',elemsL,'faces',facesL);
matL = Materials("materials.xml");

elemsR = Elements(meshR,2);
facesR = Faces(meshR);
gridR = struct('topology',meshR,'cells',elemsR,'faces',facesR);
matR = Materials("materials.xml");

% Create and set the print utility
printUtilsL = OutState(meshL,"folderName","OUT/left","timeList",[0,1,2,3,4,5,6,7,8,9,10],...
                       "writeVtk",1,"flagMatFile",1,"matFileName","OUT/left");
printUtilsR = OutState(meshR,"folderName","OUT/right","timeList",[0,1,2,3,4,5,6,7,8,9,10],...
                       "writeVtk",1,"flagMatFile",1,"matFileName","OUT/right");
% Create an object of the "Boundaries" class 
bcL = Boundaries("bcLeft.xml",gridL);
bcR = Boundaries("bcRight.xml",gridR);

% Create object handling construction of Jacobian and rhs of the model
domainL = Discretizer('Boundaries',bcL,...
                     'OutState',printUtilsL,...
                     'Materials',matL,...
                     'Grid',gridL);

domainR = Discretizer('Boundaries',bcR,...
                     'OutState',printUtilsR,...
                     'Materials',matR,...
                     'Grid',gridR);

domainL.addPhysicsSolver('StickSlipOpen.xml');
domainR.addPhysicsSolver('StickSlipOpen.xml');

domains = [domainL; domainR];

% build the mortar interface with xml input
interfaces = buildInterfaces('StickSlipOpen.xml',domains);
%interfaces = buildInterfaces('Stick.xml',domains);

solv = ActiveSetContactSolver(simParam,domains,interfaces,10);
%solv = MultidomainFCSolver(simParam,domains,interfaces);
solv.NonLinearLoop();
solv.finalizeOutput();