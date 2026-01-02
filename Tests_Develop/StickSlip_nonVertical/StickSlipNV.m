clc
close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir)

stab = "new";

% set mesh 
X = 5; Y = 10; Z = 15;
nx1 = 2; ny1 = 6; nz1 = 6;
b1 = BlockStructuredMesh([0,0.5*X;0 Y;0 Z],[nx1,ny1,nz1],1);
meshL = processGeometry(b1);

nx2 = 2; ny2 = 12; nz2 = 12;
b2 = BlockStructuredMesh([0.5*X,X;0 Y;0 Z],[nx2, ny2, nz2],1);
meshR = processGeometry(b2);

assert(mod(ny1,2) == 0 && mod(ny2,2)==0,"Number of elements along y axis " + ...
  "must be even to correctly apply symmetric bcs")

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


if stab == "old"
  fl = "_old";
elseif stab == "new"
  fl = "_new";
end



printUtilsL = OutState(meshL,"folderName",strcat("OUT/left",fl),"timeList",0:0.1:1,...
                       "writeVtk",1,"flagMatFile",1,"matFileName",strcat("OUT/left",fl));
printUtilsR = OutState(meshR,"folderName",strcat("OUT/right",fl),"timeList",0:0.1:1,...
                       "writeVtk",1,"flagMatFile",1,"matFileName",strcat("OUT/right",fl));
% Create an object of the "Boundaries" class 
setBC(Y,meshL,meshR)

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

tIni = -1;
interfaces{1}.state.traction(1:3:end) = tIni;
interfaces{1}.state.iniTraction(1:3:end) = tIni;
interfaces{1}.stateOld.iniTraction(1:3:end) = tIni;
interfaces{1}.stateOld.traction(1:3:end) = tIni;

if stab == "new"
  interfaces{1}.oldStab = false;
else
  interfaces{1}.oldStab = true;
end

interfaces{1}.outstate.VTK.setVTKFolder(strcat("OUT/Crack",fl));
interfaces{1}.outstate.matFileName = strcat("OUT/Crack",fl);


solv = ActiveSetContactSolver(simParam,domains,interfaces,10);
%solv = MultidomainFCSolver(simParam,domains,interfaces);
diary off

fname = ['log_' char(stab) '.txt'];

if exist(fname,'file')
    delete(fname)
end

diary(fname)
solv.NonLinearLoop();
diary off
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