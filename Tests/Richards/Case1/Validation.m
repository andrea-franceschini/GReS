% This file is a comparison of two approaches to describing two variables,
% the saturation and the relative permeability. The first approach proposes
% an analytical expression for these variables, while the second approach
% relies on interpolated values from a table.

close all;
% clear;
output_dir = 'Outputs';
input_dir = 'Inputs';
figures_dir = fullfile(output_dir,"Images");


%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create the Mesh object
topology = Mesh();

% Import mesh data into the Mesh object
topology.importMesh(fullfile(input_dir,'Mesh',"Column1x1x30.msh"));

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Create an object of the "Faces" class and process the face properties
faces = Faces(topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Creating boundaries conditions.
bound = Boundaries(fullfile(input_dir,'boundaries.xml'),grid);

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility
printUtils = OutState(topology,fullfile(input_dir,'outputValidation.xml'));

sol = struct();
ComparisonMaterials = ["matTable.xml","matTabular.xml"];
for sim=1:length(ComparisonMaterials)
  % Create an object of the Materials class and read the materials file
  mat = Materials(fullfile(input_dir,"Materials",ComparisonMaterials(sim)));

  % Create object handling construction of Jacobian and rhs of the model
  domain = Discretizer('Grid',grid,...
                       'Materials',mat,...
                       'Boundaries',bound,...
                       'OutState',printUtils);

  domain.addPhysicsSolver(fullfile(input_dir,'solver.xml'));

  % set initial conditions directly modifying the state object
  z = elems.mesh.cellCentroid(:,3);
  gamma_w = getFluid(mat).getFluidSpecWeight();
  wLev = 9.; % level of the water table
  domain.state.data.pressure = gamma_w*(wLev-z);

  % The modular structure of the discretizer class allow the user to easily
  % customize the solution scheme.
  % Here, a built-in fully implict solution scheme is adopted with class
  % FCSolver. This could be simply be replaced by a user defined function
  Solver = FCSolver(simParam,domain);
  % Solver = FCSolver(domain,'SaveRelError',true,'SaveBStepInf',true);

  % Solve the problem
  [simState] = Solver.NonLinearLoop();

  % Finalize the print utility
  printUtils.finalize()

  % elem vector containing elements centroid along vertical axis
  numb = 0.;
  tol = 0.01;
  nodesP = find(abs(topology.cellCentroid(:,1)-numb) < tol & abs(topology.cellCentroid(:,2)-numb) < tol);
  [~,ind] = sort(topology.cellCentroid(nodesP,3));
  nodesP = nodesP(ind);

  nrep = length(printUtils.results);
  nvars = length(printUtils.results(2).pressure);
  press = zeros(nvars,nrep);
  sw = zeros(nvars,nrep);
  t = zeros(1,nrep);
  for i=1:nrep
    press(:,i) = printUtils.results(i).pressure;
    sw(:,i) = printUtils.results(i).saturation;
    t(i) = printUtils.results(i).time;
  end
  tind = 1:length(t);

  % Getting pressure and saturation solution for specified time
  pressplot = press(nodesP,tind');

  % Vertical position of the column
  ptsZ = elems.mesh.cellCentroid(nodesP,3);
  
  % Values for normalized plots
  pos = find(ptsZ == max(ptsZ));
  H = max(topology.coordinates(:,3));

  sol(sim).pressure = pressplot./(pressplot(pos,:));
  sol(sim).position = ptsZ/H;
  sol(sim).saturation = sw(nodesP,tind);
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