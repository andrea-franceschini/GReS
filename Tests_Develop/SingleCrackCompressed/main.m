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
szBot = 3; szTop = 5; szBotF = 0.03; szTopF = 0.05; nz = 2;

com = "python Mesh/grid.py " + outFile + " " + num2str(szBot) + " " + num2str(szBotF) + " "  +...
  num2str(szTop) + " " + num2str(szTopF) + " " + num2str(nz);

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
%domains(1).grid.topology = setFaultSurface(domains(1).grid.topology);



[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);



% set verbosity 
domains(1).simparams.setVerbosity(2);

%solver = MultidomainFCSolver(domains,interfaces);

solver = ActiveSetContactSolver(domains,interfaces,10);


solver.NonLinearLoop();
solver.finalizeOutput();

%% POST PROCESSING

%sort cells based on x-coordinate
nS = solver.interfaces{2}.mesh.msh(2).nSurfaces;
cId = 2:2:nS-1;
gt = solver.interfaces{2}.slip.curr(cId);
tn = solver.interfaces{2}.traction.curr(3*cId-2);

xCoord = solver.interfaces{2}.mesh.msh(2).surfaceCentroid(cId,1)/cos(deg2rad(20));
xAnal = linspace(-1,1,1000);


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
gt_anal = K*sqrt(b^2-(b-xi).^2);

figure(1)
plot(xCoord, gt, 'k-o', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
hold on
plot(xAnal, gt_anal, 'b-', 'MarkerSize', 1, 'LineWidth', 1.5);
xlim([-1 1])
ylim([0 4.1e-3])
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\|\mathbf{g_T}\|$', 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','FontSize',14)   % <-- axis numbers in LaTeX

exportgraphics(gcf, 'gt_plot.pdf')

figure(2)
plot(xCoord, -tn, 'k-o', 'MarkerSize', 4,'MarkerFaceColor', 'k');
hold on
plot([-1 1], [tn_anal tn_anal], 'r-', 'LineWidth', 1.5)
xlim([-1 1])
ylim([-10 20])
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\sigma_n$', 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','FontSize',14)   % <-- axis numbers in LaTeX

exportgraphics(gcf, 'tn_plot.pdf')




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

