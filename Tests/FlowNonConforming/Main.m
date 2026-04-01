% clear
close all
clc


% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');
% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

% Flow non conforming model - this model does not use any additional input
% files


% Domain 1 - Top layer
mshTop = Mesh();
mshTop.importMesh('topLayer.vtk');
gridTop = struct('topology',mshTop,'cells',Elements(mshTop));
matTop = Materials();
matTop.addSolid('name',"topMat",'cellTags',1);
matTop.addFluid('specificWeight',0.0,'compressibility',4.59e-7,'dynamicViscosity',1e-6);
matTop.addPorousRock("topMat",'porosity',0.375,'permeability',5e-10,'specificWeight',2.1e1);

bcTop = Boundaries(gridTop);
bcTop.addBC('name',"topCornerPressure",...
  'type',"dirichlet",...
  'field',"node",...
  'entityListType',"bcList",...
  'entityList',8,...
  'variable',"pressure")
bcTop.addBCEvent("topCornerPressure",'time',0.0,'value',0.0);

domainTop = Discretizer('grid',gridTop,'materials',matTop,'boundaries',bcTop);

domainTop.addPhysicsSolver('SinglePhaseFlowFEM');


% Domain 2 - Bottom layer
mshBottom = Mesh();
mshBottom.importMesh('bottomLayer.vtk');
gridBottom = struct('topology',mshBottom,'cells',Elements(mshBottom));
matBottom = Materials();
matBottom.addSolid('name',"botMat",'cellTags',1);
matBottom.addFluid('specificWeight',0.0,'compressibility',4.59e-7,'dynamicViscosity',1e-6);
matBottom.addPorousRock("botMat",'porosity',0.3,'permeability',1e-10,'specificWeight',2.1e1);

bcBottom = Boundaries(gridBottom);
bcBottom.addBC('name',"bottomCornerPressure",...
  'type',"dirichlet",...
  'field',"node",...
  'entityListType',"bcList",...
  'entityList',2,...
  'variable',"pressure")
bcBottom.addBCEvent("bottomCornerPressure",'time',0.0,'value',1.0e2);

domainBottom = Discretizer('grid',gridBottom,'materials',matBottom,'boundaries',bcBottom);

domainBottom.addPhysicsSolver('SinglePhaseFlowFEM');



% Interface
domains = [domainTop,domainBottom];

interf.masterDomain = 1;
interf.masterSurface = 1;
interf.slaveDomain = 2;
interf.slaveSurface = 1;
interf.multiplierType = "dual";
interf.Quadrature.type = "SegmentBasedQuadrature";
interf.Quadrature.nGP = 4;
interfaces = InterfaceSolver.add('MeshTying',domains,interf);


simparams = SimulationParameters('simParams.xml');
printUtils = OutState('printTimes',[0.5,1,2,5,10,50],'outputFile',"flowNonConformingOut");

solver = NonLinearImplicit('simulationparameters',simparams,...
                           'domains',domains,...
                           'interface',interfaces, ...
                           'output',printUtils);
solver.simulationLoop();
