clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fname = 'constantSliding.xml';

simparams = SimulationParameters(fname);


simParam = SimulationParameters(fname);

elemsL = Elements(meshL,2);
facesL = Faces(meshL);
gridL = struct('topology',meshL,'cells',elemsL,'faces',facesL);
matL = Materials("materials.xml");

elemsR = Elements(meshR,2);
facesR = Faces(meshR);
gridR = struct('topology',meshR,'cells',elemsR,'faces',facesR);
matR = Materials("materials.xml");

% Create and set the print utility



printUtilsL = OutState(meshL,"folderName",strcat("OUT/leftBlock"),"timeList",0:20,...
                       "writeVtk",1,"flagMatFile",1,"matFileName",strcat("OUT/leftBlock"));
printUtilsR = OutState(meshR,"folderName",strcat("OUT/rightBlock"),"timeList",0:20,...
                       "writeVtk",1,"flagMatFile",1,"matFileName",strcat("OUT/rightBlock"));
% Create an object of the "Boundaries" class 
setBC(Y,meshL,meshR)

bc = Boundaries("bcLeft.xml",gridL);
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bc,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);




solver = GeneralSolver(simparams,domains,interfaces);

solver.NonLinearLoop();
solver.finalizeOutput();

% get tangential gap
gt = interfaces{1}.state.tangentialGap(1:2:end);
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

