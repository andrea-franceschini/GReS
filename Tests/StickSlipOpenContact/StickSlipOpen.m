clc
close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir)

params = readInput('StickSlipOpen.xml');

% set mesh 
X = 5; Y = 10; Z = 15;
nx1 = 2; ny1 = 6; nz1 = 6;
meshL = structuredMesh(nx1,ny1,nz1,[0,0.5*X],[0 Y],[0 Z]);

nx2 = 2; ny2 = 8; nz2 = 8;
meshR = structuredMesh(nx2,ny2,nz2,[0.5*X,X],[0 Y],[0 Z]);

assert(mod(ny1,2) == 0 && mod(ny2,2)==0,"Number of elements along y axis " + ...
  "must be even to correctly apply symmetric bcs")

% define model 

simParam = SimulationParameters(params.SimulationParameters);

elemsL = Elements(meshL,2);
facesL = Faces(meshL);
gridL = struct('topology',meshL,'cells',elemsL,'faces',facesL);
matL = Materials("materials.xml");

elemsR = Elements(meshR,2);
facesR = Faces(meshR);
gridR = struct('topology',meshR,'cells',elemsR,'faces',facesR);
matR = Materials("materials.xml");


% Create an object of the "Boundaries" class 
[bcL,bcR] = setBC(Y,gridL,gridR);




% Create object handling construction of Jacobian and rhs of the model
domainL = Discretizer('Boundaries',bcL,...
                     'Materials',matL,...
                     'Grid',gridL);

domainR = Discretizer('Boundaries',bcR,...
                     'Materials',matR,...
                     'Grid',gridR);

domainL.addPhysicsSolver('Poromechanics');
domainR.addPhysicsSolver('Poromechanics');

domains = [domainL; domainR];

% build the mortar interface with xml input
interfaces = InterfaceSolver.addInterfaces(domains,params.Interface);


% apply initial traction to the interface
tIni = -1;
interfaces{1}.state.traction(1:3:end) = tIni;
interfaces{1}.state.iniTraction(1:3:end) = tIni;
interfaces{1}.stateOld.iniTraction(1:3:end) = tIni;
interfaces{1}.stateOld.traction(1:3:end) = tIni;

printUtils = OutState("outputFike","Ouputy/StickSlipOpen","printTimes",0:20,...
                      "matFileName","Output/StickSlipOpenHistory");

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domains,...
                           'interface',interfaces, ...
                           'output',printUtils);
solver.simulationLoop();



function [bcLeft,bcRigth] = setBC(Y,gridL,gridR)

meshL = gridL.topology;
meshR = gridR.topology;

% write bc files to apply y contraint on the right node location
targetCoord = 0.5*Y;
nL = all([abs(meshL.coordinates(:,2)-targetCoord) < 1e-4,...
          abs(meshL.coordinates(:,3)) < 1e-4],2);

nR = all([abs(meshR.coordinates(:,2)-targetCoord) < 1e-4,...
          abs(meshR.coordinates(:,3)) < 1e-4],2);


bcLeft = Boundaries(gridL);
bcLeft.addBC('name',"fixBack",...
          'type',"dirichlet",...
          'field',"surface",...
          'variable',"displacements",...
          'entityListType',"tag", ...
          'entityList',5,...
          'components',"x");
bcLeft.addBCEvent("fixBack",'time',0.0,'value',0.0);

bcLeft.addBC('name',"y_bottom",...
          'type',"dirichlet",...
          'field',"node",...
          'variable',"displacements",...
          'entityListType',"bcList", ...
          'entityList',nL,...
          'components',"y");
bcLeft.addBCEvent("y_bottom",'time',0.0,'value',0.0);

bcLeft.addBC('name',"z_bottom",...
          'type',"dirichlet",...
          'field',"surface",...
          'variable',"displacements",...
          'entityListType',"tag", ...
          'entityList',1,...
          'components',"z");
bcLeft.addBCEvent("z_bottom",'time',0.0,'value',0.0);


bcRigth = Boundaries(gridR);
bcRigth.addBC('name',"x_load",...
          'type',"neumann",...
          'field',"surface",...
          'variable',"displacements",...
          'entityListType',"tag", ...
          'entityList',6,...
          'components',"x");
bcRigth.addBCEvent("x_load",'time',0.0,'value',0.0);
bcRigth.addBCEvent("x_load",'time',1.0,'value',-5.0);
bcRigth.addBCEvent("x_load",'time',6.0,'value',-5.0);
bcRigth.addBCEvent("x_load",'time',11.0,'value',0.0);
bcRigth.addBCEvent("x_load",'time',16.0,'value',0.0);
bcRigth.addBCEvent("x_load",'time',20.0,'value',1.0);

bcRigth.addBC('name',"y_bottom",...
          'type',"dirichlet",...
          'field',"node",...
          'variable',"displacements",...
          'entityListType',"bcList", ...
          'entityList',nR,...
          'components',"y");
bcRigth.addBCEvent("y_bottom",'time',0.0,'value',0.0);

bcRigth.addBC('name',"z_bottom",...
          'type',"dirichlet",...
          'field',"surface",...
          'variable',"displacements",...
          'entityListType',"tag", ...
          'entityList',1,...
          'components',"z");
bcRigth.addBCEvent("z_bottom",'time',0.0,'value',0.0);




end