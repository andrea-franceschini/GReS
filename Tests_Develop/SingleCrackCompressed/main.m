clear
close all
clc

 
scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);


fprintf('Single crack under compression \n')
fprintf('___________________\n\n')

% write mesh files
outFile = "SingleCrack_grid";
NXb = 40; NYb = 10; NXt = 80; NYt = 7; nz = 4;

com = "python Mesh/domain.py " + outFile + " " + num2str(NXb) + " " + num2str(NYb) + " "  +...
  num2str(NXt) + " " + num2str(NYt) + " " + num2str(nz);

system(com);

grid = Mesh();
grid.importMesh(outFile + ".vtk");
% 

plotFunction(grid,'singleCrack_grid',zeros(grid.nNodes,1));

% write BC files
setBCfiles(grid);

domainFile = 'Domains/domains.xml';
interfFile = 'Domains/interface.xml';

domains = buildModel(domainFile); 


% define the internal fault to be approx m long
domains(1).grid.topology = setFaultSurface(domains(1).grid.topology);



[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);



% set verbosity 
domains(1).simparams.setVerbosity(2);

%solver = MultidomainFCSolver(domains,interfaces);

solver = ActiveSetContactSolver(domains,interfaces,10);


solver.NonLinearLoop();
solver.finalizeOutput();




function msh = setFaultSurface(msh)
  % get nodes belonging to current surface 3
  surfInFault = find(msh.surfaceTag==9);
  n = unique(msh.surfaces(surfInFault,:));
  xTarget = 1;
  nInFault = abs(msh.coordinates(n,1))<xTarget;
  nInFault = n(nInFault);
  isNotSurfInFault = ~all(ismember(msh.surfaces(surfInFault,:),nInFault),2);
  isNotSurfInFault = surfInFault(isNotSurfInFault);
  msh.surfaceTag(isNotSurfInFault) = msh.nSurfaceTag + 1;
  msh.nSurfaceTag = msh.nSurfaceTag + 1; 
end

%% plot profiles of multipliers along vertical axis (avoid opening paraview)
%plotStep(solver.results,2);






% 

