clc
close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);


cd(scriptDir)

fname = "StickBubble.xml";


% set mesh 
b1 = BlockStructuredMesh([0,2.5;0 10;0 15],[3 8 8],1);
meshL = processGeometry(b1);

b2 = BlockStructuredMesh([2.5,5;0 10;0 15],[3 16 16],1);
meshR = processGeometry(b2);

% define model 

simParam = SimulationParameters(fname);

elemsL = Elements(meshL,3);
facesL = Faces(meshL);
gridL = struct('topology',meshL,'cells',elemsL,'faces',facesL);
matL = Materials("materials.xml");

elemsR = Elements(meshR,3);
facesR = Faces(meshR);
gridR = struct('topology',meshR,'cells',elemsR,'faces',facesR);
matR = Materials("materials.xml");

% Create and set the print utility
printUtilsL = OutState(meshL,"folderName","OUT/left","timeList",[0,1,2,3,4,5,6,7,8,9,10],...
                       "writeVtk",1,"flagMatFile",1,"matFileName","OUT/left");
printUtilsR = OutState(meshR,"folderName","OUT/right","timeList",[0,1,2,3,4,5,6,7,8,9,10],...
                       "writeVtk",1,"flagMatFile",1,"matFileName","OUT/right");

% Create an object of the "Boundaries" class 
setBC(10,meshL,meshR);
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

domainL.addPhysicsSolver(fname);
domainR.addPhysicsSolver(fname);

domains = [domainL; domainR];

% build the mortar interface with xml input
interfaces = buildInterfaces(fname,domains);
%interfaces = buildInterfaces('Stick.xml',domains);

tIni = -1;
interfaces{1}.state.multipliers(1:3:end) = tIni;
interfaces{1}.state.iniMultipliers(1:3:end) = tIni;
interfaces{1}.stateOld.multipliers(1:3:end) = tIni;
interfaces{1}.stateOld.iniMultipliers(1:3:end) = tIni;


%solv = ActiveSetContactSolver(simParam,domains,interfaces,10);
solv = MultidomainFCSolver(simParam,domains,interfaces);
solv.NonLinearLoop();
solv.finalizeOutput();


function setBC(Y,meshL,meshR)

% write bc files to apply y contraint on the right node location
targetCoord = 0.5*Y;
nL = all([abs(meshL.coordinates(:,2)-targetCoord) < 1e-4,...
          abs(meshL.coordinates(:,3)) < 1e-4],2);

nR = all([abs(meshR.coordinates(:,2)-targetCoord) < 1e-4,...
          abs(meshR.coordinates(:,3)) < 1e-4],2);

strBCL = readstruct("bcLeft.xml",AttributeSuffix = "");
strBCL.BC(2).BCentities.bcList = find(nL);
writestruct(strBCL,"bcLeft.xml",AttributeSuffix="");

strBCR = readstruct("bcRight.xml",AttributeSuffix = "");
strBCR.BC(2).BCentities.bcList = find(nR);
writestruct(strBCR,"bcRight.xml",AttributeSuffix="");


end