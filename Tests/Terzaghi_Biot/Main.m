close all;
% clear;
output_dir = 'Output';
input_dir = 'Input';
figures_dir = fullfile(output_dir,"Images");

% input_dir = 'Input';
% figures_dir = 'Images';

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
% cd(scriptDir);

%% ------------------------------------------------------------------------

% Set parameters of the simulation
simParam = SimulationParameters(fullfile(scriptDir,input_dir,"simparam.xml"));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(scriptDir,input_dir,"materials.xml"));

%% ------------------------------ Set up the Domain -----------------------
% Create the Mesh object
topology = Mesh();

% Import mesh data into the Mesh object
MeshFile = ["Column.msh","Column_hexa.msh","Column_tetra.msh"];
topology.importMesh(fullfile(scriptDir,input_dir,"Mesh",MeshFile(2)));

% Create an object of the "Elements" class and process the element properties
gaussOrder = 2;
elems = Elements(topology,gaussOrder);

% Create an object of the "Faces" class and process the face properties
faces = Faces(topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);


% Creating boundaries conditions.
bound = Boundaries(fullfile(scriptDir,input_dir,"boundaries.xml"),grid);

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility for the solution
printUtils = OutState(topology,fullfile(scriptDir,input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

domain.addPhysicsSolver('solver_TPFA.xml');

% In this version of the code, the user can assign initial conditions only
% manually, by directly modifying the entries of the state structure. 
% In this example, we use a user defined function to apply Terzaghi initial
% conditions to the state structure
F = -10; % vertical force
applyTerzaghiIC(domain.state,mat,topology,F);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = GeneralSolver(simParam,domain);

% Solve the problem
Solver.NonLinearLoop();

% Finalize the print utility
domain.outstate.finalize()

% calling analytical solution script
Terzaghi_analytical(topology, mat, abs(F),[15,30,60,90,120,180],output_dir)


%% --------------------- Post Processing the Results ----------------------
if true
  image_dir = fullfile(pwd,figures_dir);
  if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(figures_dir)
  else
    mkdir(figures_dir)
  end

  load(fullfile(output_dir,"Terzaghi_Analytical.mat"))
  %Post processing using MAT-FILE
  % obtain vector contain list of nodes along vertical axis (with x,y=0)
  nodesU = find(topology.coordinates(:,1)+topology.coordinates(:,2)==0);
  [~,ind] = sort(topology.coordinates(nodesU,3));
  nodesU = nodesU(ind);


  % elem vector containing elements centroid along vertical axis
  flowscheme = getPhysicsSolver(domain,"BiotFullySaturated").getFlowScheme();
  if strcmp(flowscheme,"FEM")
    nodesP = nodesU;
  else
    nodesP = find(topology.cellCentroid(:,1) + topology.cellCentroid(:,2) < 0.51);
    [~,ind] = sort(topology.cellCentroid(nodesP,3));
    nodesP = nodesP(ind);
  end

  %Getting pressure and displacement solution for specified time from MatFILE
  press = [printUtils.results.pressure];
  disp = [printUtils.results.displacements];
  pressplot = press(nodesP,1:end);
  dispplot = disp(3*nodesU,1:end);

  H = max(topology.coordinates(:,3));
  p0 = max(press(:,1));

  %Plotting solution
  if strcmp(flowscheme,"FVTPFA")
    ptsY = topology.cellCentroid(nodesP,3);
  else
    ptsY = topology.coordinates(nodesP,3);
  end

  figure('Position', [100, 100, 700, 700])
  hold on
  plotObj1 = plot(pressplot/p0,ptsY/H,'k.', 'LineWidth', 1, 'MarkerSize', 15);  
  plotObj2 = plot(p/p0,z/H,'k-', 'LineWidth', 1);
  xlabel('p/p_0')
  ylabel('z/H')
  legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'northeast');
  %title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')
  axis tight
  xlim([0 1.02])

  set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10, 'XGrid', 'on', 'YGrid', 'on')

  % export figure with quality
  stmp = fullfile(figures_dir,'Terzaghi_pressure.png');
  exportgraphics(gcf,stmp,'Resolution',400)

  figure('Position', [100, 100, 700, 700])
  hold on
  plotObj1 = plot(-dispplot/H,topology.coordinates(nodesU,3)/H,'k.', 'LineWidth', 1, 'MarkerSize', 15);
  plotObj2 = plot(u/H,z/H,'k-',  'LineWidth', 1);
  xlabel('u_z/H')
  ylabel('z/H')
  %title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')
  legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'southeast');

  set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 10, 'XGrid', 'on', 'YGrid', 'on')

  % export figure with quality
  stmp = fullfile(figures_dir,'Terzaghi_disp.png');
  exportgraphics(gcf,stmp,'Resolution',400)
end