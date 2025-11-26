%% write mesh files
outFile = "SingleCrack_grid";

grid = Mesh();
grid.importMesh(outFile + ".vtk");
% 

domainFile = 'Domains/domains.xml';
interfFile = 'Domains/interface.xml';

domains = buildModel(domainFile); 

% set verbosity 
domains(1).simparams.setVerbosity(0);

[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);

solver = ActiveSetContactSolver(domains,interfaces,10);

solver.NonLinearLoop();
solver.finalizeOutput();

%% Validate 
nS = solver.interfaces{2}.mesh.msh(2).nSurfaces;
cId = 2:2:nS-1;
gt = solver.interfaces{2}.slip.curr(cId);
tn = solver.interfaces{2}.traction.curr(3*cId-2);

xCoord = solver.interfaces{2}.mesh.msh(2).surfaceCentroid(cId,1)/cos(deg2rad(20));

% analytical solutions
b = 1;
% real angle of the generated fault
c = solver.interfaces{2}.mesh.msh(2).surfaceCentroid(end,:);
psi = atan(abs(c(2)/c(1)));
sigma = 100;
nu = 0.25;
E = 15000;
theta = deg2rad(30);
tn_anal = sigma*(sin(psi))^2;
xi = linspace(0,2,1000);
K = 4*(1-nu^2)*(sigma*sin(psi)*(cos(psi)-sin(psi)*tan(theta)))/E;
gt_anal = K*sqrt(b^2-(b-(xCoord+1)).^2);
gt_anal = flip(gt_anal);

% compare normal traction
ntn = numel(tn);
ntn_del = round(ntn/10);
err_tn = norm(tn(ntn_del:end-ntn_del)+tn_anal);
assert(err_tn < 1e0,"Normal traction not validated")

% compare tangential gap
err_gt = norm(gt-gt_anal);
assert(err_gt < 1e-2,"Tangential gap not validated")

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







% 

