close all;
% clear;
input_dir = 'Input/';
figures_dir = 'Output/Figs/';

typeDiscretization = "FVTPFA";
solverName = strcat("SinglePhaseFlow",typeDiscretization);

% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% output parameters
printUtils = OutState('outputFile',"Output/results",'printTimes',[1,2,3,4,5,6,7,8,9,10]);

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create the Mesh object
grid = Grid();

% Import mesh data into the Mesh object
grid.importMesh(fullfile(input_dir,'domain.msh'));

% Creating boundaries conditions.
bound = Boundaries(grid,fullfile(input_dir,'boundaries.xml'));


% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound);

domain.addPhysicsSolver(solverName,'steadyState',0);


% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = 1.e5;

solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

%% --------------------- Post Processing the Results ----------------------
postproc=true;
if postproc
  image_dir = fullfile(pwd,figures_dir);
  if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(figures_dir)
  else
    mkdir(figures_dir)
  end

  % Saving a temporary variabel.
  pressure = [printUtils.results.pressure];

  % Ajusting the time position.
  t = [printUtils.results.time];
  tind = 2:length(t);
  t_max = t(end);
  t = t(tind);
  tstr = strcat(num2str(t'),' second');

  center = grid.cells.center;

  %Getting pressure and saturation solution for specified time from MatFILE
  numbA = 5.;
  numbB = 5.;
  switch typeDiscretization
    case "FEM"
      tol = 1e-3;
      nodesP = find(abs(grid.coordinates(:,1)-numbA) < tol & abs(grid.coordinates(:,3)-numbB) < tol);
      [~,ind] = sort(grid.coordinates(nodesP,2));
      nodesP = nodesP(ind);

      pressplot = pressure(nodesP);

      % Values for normalized plots
      H = max(grid.coordinates(nodesP,2));

      % Location a column to be the plot position.
      pts = grid.coordinates(nodesP,2)/H;
    case "FVTPFA"
      tol = 0.4;
      nodesP = find(abs(center(:,1)-numbA) < tol & abs(center(:,3)-numbB) < tol);
      [~,ind] = sort(center(nodesP,2));
      nodesP = nodesP(ind);

      pressplot = pressure(nodesP);

      % Values for normalized plots
      H = max(center(nodesP,2));

      % Location a column to be the plot position.
      pts = center(nodesP,2)/H;
    otherwise
  end

  %Plotting pressure head
  figure(1)
  plot(pts,pressplot,'.-', 'LineWidth', 2, 'MarkerSize', 14);
  hold on
  ylabel('pressure (Pa)')
  xlabel('distance (%)')
  legend(tstr)
  set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
  % export figure with quality
  stmp = fullfile(image_dir,'pressure.png');
  exportgraphics(gcf,stmp,'Resolution',400)
end