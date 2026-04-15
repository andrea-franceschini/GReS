clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fname = fullfile('Input','pressurizedCrack.xml');

params = readInput(fname);

simparams = SimulationParameters(params.SimulationParameters);

% coordinate arrays
% Domain sizes
X = 40.0;
Y = 1;
Z = X;

grid = structuredMesh([10,201,10],1,[10,61,10],[-0.5*X,-2,2,0.5*X],[-0.5*Y 0.5*Y],[-0.5*Z,-2,2,0.5*Z]);

%assert(3*mesh.nNodes < 2e5,"Mesh is too fine")

%%

mat = Materials(params.Materials);

printUtils = OutState("outputFile",strcat("Output/Sneddon"),"printTimes",1,...
                      "matFileName",strcat("Output/Sneddon"));

bc = Boundaries(grid);

bc.addBC('name',"z_fix",...
          'type',"dirichlet",...
          'field',"surface",...
          'variable',"displacements",...
          'entityListType',"tag", ...
          'entityList',[1,2],...
          'components',"z");
bc.addBCEvent("z_fix",'time',0.0,'value',0.0);

bc.addBC('name',"y_fix",...
          'type',"dirichlet",...
          'field',"surface",...
          'variable',"displacements",...
          'entityListType',"tag", ...
          'entityList',[3,4],...
          'components',"x");
bc.addBCEvent("y_fix",'time',0.0,'value',0.0);

bc.addBC('name',"x_fix",...
          'type',"dirichlet",...
          'field',"surface",...
          'variable',"displacements",...
          'entityListType',"tag", ...
          'entityList',[5,6],...
          'components',"x");
bc.addBCEvent("x_fix",'time',0.0,'value',0.0);



% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bc,...
                     'Materials',mat,...
                     'Grid',grid);


domain.addPhysicsSolvers(params.Solver);

% set neumann traction in embedded fracture
efem = getPhysicsSolver(domain,"EmbeddedFractureMechanics");
efem.bcTraction(1:3:end) = -2;

solver = NonLinearImplicit('simulationparameters',simparams,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();


%% analytical solution processing


% fracture angle
efem = getPhysicsSolver(domain,"EmbeddedFractureMechanics");
f = efem.fractureMesh.surfaces;
% direction along fracture
v = f.tang1(1,:);
psi = acos(v*[1;0;0]);
angle = rad2deg(psi);

fractureSize = 2.0;

% get unique set of cut cells along the fracture
[v,i] = sort(f.center(:,2),"ascend");
id = abs(diff(v)) > 1e-3;
centers = f.center(i(1:find(id)),:);

P1 = [-fractureSize/2,0,0];

% lenght of the fault
L = fractureSize/cos(deg2rad(angle));
xi = f.center(:,1)/cos(deg2rad(angle));


gn = efem.domain.state.data.fractureJump(1:3:end); 


% analytical solutions
b = L/2;

sigma = 2;

nu = 0.25;
E = 25000;
theta = deg2rad(30);
% tn_anal = sigma*(sin(psi))^2;
xi_anal = linspace(-L/2,L/2,1000);
K = 4*(1-nu^2)*sigma/E;
%K = 4*(1-nu^2)*(sigma*sin(psi)*(cos(psi)-sin(psi)*tan(theta)))/E;
gn_anal = K*sqrt(b^2-xi.^2);
% gt_anal = flip(gt_anal);

err = norm(1e3*f.area(1)*(gn-gn_anal));
fprintf("Error norm: %1.4e \n",err)

% plot

% gn_anal_plot = K*sqrt(b^2-xi_anal.^2);
% figure(1)
% plot(xi, gn, 'k-o', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
% hold on
% plot(xi_anal, gn_anal_plot, 'b-', 'MarkerSize', 1, 'LineWidth', 1.5);
% xlim([-1.5*b 1.5*b])
% ylim([0 1e-3])
% xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\|\mathbf{g_T}\|$', 'Interpreter', 'latex', 'FontSize', 14)
% set(gca,'TickLabelInterpreter','latex','FontSize',14)   % <-- axis numbers in LaTeX

%exportgraphics(gcf, fullfile('Output','gn_plot.png'))


% 

