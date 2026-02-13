clear
close all

scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

cd(scriptDir);

fname = 'pressurizedCrack.xml';

simparams = SimulationParameters(fname);

% coordinate arrays
% Domain sizes
X = 40.0;
Y = 1;
Z = X;

mesh = structuredMesh(521,1,101,[-0.5*X 0.5*X],[-0.5*Y 0.5*Y],[-0.5*Z 0.5*Z]);

%assert(3*mesh.nNodes < 2e5,"Mesh is too fine")

%%
elems = Elements(mesh,2);
faces = Faces(mesh);
grid = struct('topology',mesh,'cells',elems,'faces',faces);
mat = Materials(fname);


printUtils = OutState("folderName",strcat("Output/Sneddon"),"timeList",1,...
                       "writeVtk",1,"flagMatFile",1,"matFileName",strcat("Output/Sneddon"));


bc = Boundaries(fname,grid);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bc,...
                     'Materials',mat,...
                     'Grid',grid);


domain.addPhysicsSolver(fname);

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

% direction along fracture
v = efem.cutTang1(1,:);
psi = acos(v*[1;0;0]);
angle = rad2deg(psi);

fractureSize = 2.0;

% get unique set of cut cells along the fracture
[v,i] = sort(efem.cutCenters(:,2),"ascend");
id = abs(diff(v)) > 1e-3;
centers = efem.cutCenters(i(1:find(id)),:);

P1 = [-fractureSize/2,0,0];

% lenght of the fault
L = fractureSize/cos(deg2rad(angle));
xi = efem.cutCenters(:,1)/cos(deg2rad(angle));


gn = efem.domain.state.data.fractureJump(1:3:end); 


% analytical solutions
b = L/2;

sigma = 2;

nu = 0.25;
E = 15000;
theta = deg2rad(30);
% tn_anal = sigma*(sin(psi))^2;
xi_anal = linspace(-L/2,L/2,1000);
K = 4*(1-nu^2)*sigma/E;
%K = 4*(1-nu^2)*(sigma*sin(psi)*(cos(psi)-sin(psi)*tan(theta)))/E;
gn_anal = K*sqrt(b^2-xi.^2);
% gt_anal = flip(gt_anal);

%%

% err_gt = norm(gt-gt_anal);
% assert(err_gt < 1e-2,"Tangential gap not validated")

%% plot

gn_anal_plot = K*sqrt(b^2-xi_anal.^2);
figure(1)
plot(xi, gn, 'k-o', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
hold on
plot(xi_anal, gn_anal_plot, 'b-', 'MarkerSize', 1, 'LineWidth', 1.5);
xlim([-1.5*b 1.5*b])
ylim([0 1e-3])
xlabel('$\xi$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\|\mathbf{g_T}\|$', 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'TickLabelInterpreter','latex','FontSize',14)   % <-- axis numbers in LaTeX

%exportgraphics(gcf, 'gn_plot.pdf')




% 

