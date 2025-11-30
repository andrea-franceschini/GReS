clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fileName = "singleCrackCompressed.xml";

%% Run model

simparams = SimulationParameters(fileName);

[domains,interfaces] = buildModel(fileName); 

solver = ActiveSetContactSolver(simparams,domains,interfaces,10);

solver.NonLinearLoop();
solver.finalizeOutput();

%% Validate
nS = getMesh(interfaces{2},MortarSide.slave).nSurfaces;
cId = 2:2:nS-1;
gt = interfaces{2}.state.slip(cId);
tn = interfaces{2}.state.traction(3*cId-2);

xCoord = getMesh(interfaces{2},MortarSide.slave).surfaceCentroid(cId,1)/cos(deg2rad(20));

% analytical solutions
b = 1;
% real angle of the generated fault
c = getMesh(interfaces{2},MortarSide.slave).surfaceCentroid(end,:);
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

%% Check normal traction
ntn = numel(tn);
ntn_del = round(ntn/10);
err_tn = norm(tn(ntn_del:end-ntn_del)+tn_anal);
assert(err_tn < 1e0,"Normal traction not validated")

%% Checl tangential gap
err_gt = norm(gt-gt_anal);
assert(err_gt < 1e-2,"Tangential gap not validated")
