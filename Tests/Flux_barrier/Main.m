close all;
% clear;
input_dir = 'Inputs/';
figures_dir = 'Figs/';

typeDiscretization = "FVTPFA";
%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create the Mesh object
topology = Mesh();

% Import mesh data into the Mesh object
topology.importMesh(fullfile(input_dir,'domain.msh'));

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Create an object of the "Faces" class and process the face properties
faces = Faces(topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Creating boundaries conditions.
bound = Boundaries(fullfile(input_dir,'boundaries.xml'),grid);

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility for the solution
printUtils = OutState(topology,fullfile(input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound,...
                     'OutState',printUtils);

switch typeDiscretization
  case "FEM"
    domain.addPhysicsSolver(fullfile(input_dir,'solverFEM.xml'));
  case "FVTPFA"
    domain.addPhysicsSolver(fullfile(input_dir,'solverFVTPFA.xml'));
end

% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = 1.e5;
% domain.state.data.potential(:) = domain.state.data.pressure+ mat.getFluid().getFluidSpecWeight()*topology.cellCentroid(:,3);

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

  %Getting pressure and saturation solution for specified time from MatFILE
  numbA = 5.;
  numbB = 5.;
  switch typeDiscretization
    case "FEM"
      tol = 1e-3;
      nodesP = find(abs(topology.coordinates(:,1)-numbA) < tol & abs(topology.coordinates(:,3)-numbB) < tol);
      [~,ind] = sort(topology.coordinates(nodesP,2));
      nodesP = nodesP(ind);

      pressplot = pressure(nodesP);

      % Values for normalized plots
      H = max(topology.coordinates(nodesP,2));

      % Location a column to be the plot position.
      pts = topology.coordinates(nodesP,2)/H;
    case "FVTPFA"
      tol = 0.4;
      nodesP = find(abs(topology.cellCentroid(:,1)-numbA) < tol & abs(topology.cellCentroid(:,3)-numbB) < tol);
      [~,ind] = sort(topology.cellCentroid(nodesP,2));
      nodesP = nodesP(ind);

      pressplot = pressure(nodesP);

      % Values for normalized plots
      H = max(topology.cellCentroid(nodesP,2));

      % Location a column to be the plot position.
      pts = topology.cellCentroid(nodesP,2)/H;
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