close all;
% clear;
input_dir = 'Inputs';
figures_dir = 'Figs';

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,"Materials",'matTable.xml'));

% Create the Mesh object
topology = Mesh();

% Choosing the mesh file
availMesh = [ "Column.msh", "Column1x1x30.msh", "Column4x4x40.msh"];
fileName = availMesh(2);

% Import mesh data into the Mesh object
topology.importMesh(fullfile(input_dir,'Mesh',fileName));

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
printUtils = OutState(topology,fullfile(input_dir,'output.xml'));

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

%% --------------------- Post Processing the Results ----------------------
if true
  image_dir = fullfile(pwd,figures_dir);
  if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(figures_dir)
  else
    mkdir(figures_dir)
  end

  % elem vector containing elements centroid along vertical axis
  numb = 0.;
  % numb = 0.125;
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
  t_max = t(end);
  t = t(tind)/t_max;
  tstr = strcat(num2str(t'),' T');

  % Getting pressure and saturation solution for specified time
  pressplot = press(nodesP,tind');
  swplot = sw(nodesP,tind);

  % Vertical position of the column
  ptsZ = elems.mesh.cellCentroid(nodesP,3);
  
  % Values for normalized plots
  pos = find(ptsZ == max(ptsZ));
  H = max(topology.coordinates(:,3));
  ptsZ = ptsZ/H;

  figure('Position', [100, 100, 700, 700])
  hold on
  plot(pressplot./(pressplot(pos,:)),ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
  xlabel('p/p_{top}')
  ylabel('z/H')
  legend(tstr, 'Location', 'northwest')
  % legend(tstr, 'Location', 'southeast')
  set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
  % export figure with quality
  stmp = fullfile(image_dir,'pressure.png');
  exportgraphics(gcf,stmp,'Resolution',400)

  figure('Position', [100, 100, 700, 700])
  plot(swplot,ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
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