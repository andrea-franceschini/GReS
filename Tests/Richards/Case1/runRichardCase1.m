close all;
% clear;
output_dir = 'Output';
input_dir = 'Input';
figures_dir = fullfile(output_dir,"Images");


%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,"Materials",'matTable.xml'));

% Create the Mesh object
grid = Grid();

% Choosing the mesh file
availMesh = [ "Column.msh", "Column1x1x30.msh", "Column4x4x40.msh"];
fileName = availMesh(2);

% Import mesh data into the Mesh object
grid.importMesh(fullfile(input_dir,'Mesh',fileName));

% Creating boundaries conditions.
bound = Boundaries(grid,fullfile(input_dir,'boundaries.xml'));

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility
printUtils = OutState(fullfile(input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound);

domain.addPhysicsSolver('VariablySaturatedFlow');

% set initial conditions directly modifying the state object
z = grid.cells.center(:,3);
gamma_w = getFluid(mat).getSpecificWeight();
wLev = 9.; % level of the water table
p = gamma_w*(wLev-z);
setState(domain,p,"pressure");

% Set and solve the simulation
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

%% --------------------- Post Processing the Results ----------------------
if true
  image_dir = fullfile(pwd,figures_dir);
  if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(figures_dir)
  else
    mkdir(figures_dir)
  end

  center = grid.cells.center;

  % elem vector containing elements centroid along vertical axis
  numb = 0.;
  % numb = 0.125;
  tol = 0.01;
  nodesP = find(abs(center(:,1)-numb) < tol & abs(center(:,2)-numb) < tol);
  [~,ind] = sort(center(nodesP,3));
  nodesP = nodesP(ind);

  tstr = strcat(num2str(printUtils.timeList'),' T');
  nrep = length(printUtils.timeList);
  nvars = length(nodesP);
  press = zeros(nvars,nrep);
  sw = zeros(nvars,nrep);
  for i=1:length(printUtils.timeList)
    press(:,i) = printUtils.results(i).pressure(nodesP);
    sw(:,i) = printUtils.results(i).saturation(nodesP);
  end

  % Vertical position of the column
  ptsZ = center(nodesP,3);
  
  % Values for normalized plots
  pos = find(ptsZ == max(ptsZ));
  H = max(grid.coordinates(:,3));
  ptsZ = ptsZ/H;

  figure('Position', [100, 100, 700, 700])
  hold on
  plot(press./(press(pos,:)),ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
  xlabel('p/p_{top}')
  ylabel('z/H')
  legend(tstr, 'Location', 'northwest')
  % legend(tstr, 'Location', 'southeast')
  set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
  % export figure with quality
  stmp = fullfile(image_dir,'pressure.png');
  exportgraphics(gcf,stmp,'Resolution',400)

  figure('Position', [100, 100, 700, 700])
  plot(sw,ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
  hold on
  xlabel('S_w')
  ylabel('z/H')
  legend(tstr, 'Location', 'northeast')
  % legend(tstr, 'Location', 'southwest')
  str = strcat('t = ',tstr);
  set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
  % export figure with quality
  stmp = fullfile(image_dir,'saturation.png');
  exportgraphics(gcf,stmp,'Resolution',400)
end