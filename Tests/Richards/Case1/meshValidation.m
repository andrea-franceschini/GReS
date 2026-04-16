% This file is a comparison of two approaches to describing two variables,
% the saturation and the relative permeability. The first approach proposes
% an analytical expression for these variables, while the second approach
% relies on interpolated values from a table.

close all;
% clear;
output_dir = 'Output';
input_dir = 'Input';
figures_dir = fullfile(output_dir,"Images");


%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create the Mesh object
grid = Grid();

% Import mesh data into the Mesh object
grid.importMesh(fullfile(input_dir,'Mesh',"Column1x1x30.msh"));

% Creating boundaries conditions.
bound = Boundaries(grid,fullfile(input_dir,'boundaries.xml'));

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility
printUtils = OutState(fullfile(input_dir,'outputValidation.xml'));

sol = struct();
ComparisonMaterials = ["matTable.xml","matTabular.xml"];


for sim=1:length(ComparisonMaterials)
  % Create an object of the Materials class and read the materials file
  mat = Materials(fullfile(input_dir,"Materials",ComparisonMaterials(sim)));

  % Create object handling construction of Jacobian and rhs of the model
  domain = Discretizer('Grid',grid,...
                       'Materials',mat,...
                       'Boundaries',bound);

  domain.addPhysicsSolver('VariablySaturatedFlow');

  % set initial conditions directly modifying the state object
  z = elems.mesh.cellCentroid(:,3);
  gamma_w = getFluid(mat).getSpecificWeight();
  wLev = 9.; % level of the water table
  domain.state.data.pressure = gamma_w*(wLev-z);

  % Set and solve the simulation
  solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);  
  solver.simulationLoop();

  % elem vector containing elements centroid along vertical axis
  numb = 0.;
  tol = 0.01;
  nodesP = find(abs(grid.cellCentroid(:,1)-numb) < tol & abs(grid.cellCentroid(:,2)-numb) < tol);
  [~,ind] = sort(grid.cellCentroid(nodesP,3));
  nodesP = nodesP(ind);

  nrep = length(printUtils.timeList);
  nvars = length(nodesP);
  press = zeros(nvars,nrep);
  sw = zeros(nvars,nrep);
  for i=1:length(printUtils.timeList)
    press(:,i) = printUtils.results(i).pressure(nodesP);
    sw(:,i) = printUtils.results(i).saturation(nodesP);
  end

  % Vertical position of the column
  ptsZ = elems.mesh.cellCentroid(nodesP,3);
  
  % Values for normalized plots
  pos = find(ptsZ == max(ptsZ));
  H = max(grid.coordinates(:,3));

  sol(sim).pressure = press./(press(pos,:));
  sol(sim).position = ptsZ/H;
  sol(sim).saturation = sw;
end
clearvars -except sol figures_dir

%% --------------------- Post Processing the Results ----------------------
image_dir = fullfile(pwd,figures_dir);
if ~isfolder(image_dir)
  mkdir(figures_dir)
end

t = [10 50 100];
t_max = t(end);
t = t/t_max;
tAnalytical = strcat("Ana - ",num2str(t',"%.1f")," T");
tTabular = strcat("Tab - ",num2str(t',"%.1f")," T");
tstr = strcat([tAnalytical; tTabular]);

figure('Position', [100, 100, 700, 700])
hold on
plot(sol(1).pressure,sol(1).position,'k-', 'LineWidth', 2, 'MarkerSize', 14);
plot(sol(2).pressure,sol(2).position,'k.', 'LineWidth', 2, 'MarkerSize', 14);
xlabel('p/p_{top}')
ylabel('z/H')
legend(tstr, 'Location', 'southeast', 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = fullfile(image_dir,'pressure_validation.png');
exportgraphics(gcf,stmp,'Resolution',400)

figure('Position', [100, 100, 700, 700])
hold on
plot(sol(1).saturation,sol(1).position,'k-', 'LineWidth', 2, 'MarkerSize', 14);
plot(sol(2).saturation,sol(2).position,'k.', 'LineWidth', 2, 'MarkerSize', 14);
xlabel('S_w')
ylabel('z/H')
legend(tstr, 'Location', 'southwest', 'NumColumns', 2)
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = fullfile(image_dir,'saturation_validation.png');
exportgraphics(gcf,stmp,'Resolution',400)