close all;
% clear;
output_dir = 'Output';
input_dir = 'Input';
figures_dir = fullfile(output_dir,"Images");

profile on

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
% cd(scriptDir);

%% ------------------------------------------------------------------------

type = "FVTPFA";
solverName = strcat("SinglePhaseFlow",type);

% Set parameters of the simulation
simParam = SimulationParameters(fullfile(scriptDir,input_dir,"simparam.xml"));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(scriptDir,input_dir,"materials.xml"));

%% ------------------------------ Set up the Domain -----------------------
% Create the Mesh object
grid = structuredMesh(20,2,20,[0,1],[0,0.1],[0,1]);
% grid = Grid();
% 
% % Set the mesh input file name
% grid.importMesh(fullfile(scriptDir,input_dir,"Mandel_Mesh.msh"));


% Creating boundaries conditions.
bound = Boundaries(grid,fullfile(scriptDir,input_dir,"boundaries.xml"));

%% ------------------ Set up and Calling the Solver -----------------------
% Create and set the print utility
printUtils = OutState('outputFile',"Output/results",'printTimes',[0.05,1,2.5,5]);

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Boundaries',bound,...
                     'Materials',mat,...
                     'Grid',grid);

solverIn.(solverName) = [];

domain.addPhysicsSolver('BiotFullyCoupled',solverIn)

% In this version of the code, the user can assign initial conditions only
% manually, by directly modifying the entries of the state structure. 
% In this example, we use a user defined function to apply Terzaghi initial
% conditions to the state structure
F = -10; % vertical force
state = applyMandelIC(domain.state,mat,grid,F);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

% calling analytical solution script
Mandel_Analytical(grid, mat, abs(F),[0.05,1,2.5,5],output_dir)

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

  %Post processing using MAT-FILE
  %list of nodes along vertical axis (with x,y=0)
  tol = 0.001;
  elemP1 = find(abs(center(:,3) - 0.025) < tol);
  elemP2 = find(abs(center(:,2) - 0.025) < tol);
  elemP = intersect(elemP1, elemP2);
  nodesX1 = find(abs(grid.coordinates(:,2)-0.05)<tol) ;
  nodesX2 = find(abs(grid.coordinates(:,3)-0.7)<tol);
  nodesX = intersect(nodesX1,nodesX2);
  nodesZ1 = find(abs(grid.coordinates(:,1)-0.6)<tol);
  nodesZ2 = find(abs(grid.coordinates(:,2)-0.05)<tol);
  nodesZ = intersect(nodesZ1,nodesZ2);
  [coordsP,ind] = sort(center(elemP,1));
  elemP = elemP(ind);
  [coordsX,ind] = sort(grid.coordinates(nodesX,1));
  nodesX = nodesX(ind);
  [coordsZ,ind] = sort(grid.coordinates(nodesZ,3));
  nodesZ = nodesZ(ind);

  %Getting pressure and displacement solution for specified output times from MatFILE
  press = [printUtils.results.pressure];
  disp = [printUtils.results.displacements];
  pressNum = press(elemP,1:end);
  dispXNum = disp(3*nodesX-2,1:end);
  dispZNum = disp(3*nodesZ,1:end);

  % Loading analytical solution MAT-FILE into workspace
  load(fullfile(output_dir,"Mandel_Analytical.mat"));

  % Plotting solution
  % Pressure
  figure('Position', [100, 100, 700, 700])
  hold on
  plotObj1 = plot(center(elemP,1),pressNum,'k.', 'LineWidth', 1, 'MarkerSize', 15);
  plotObj2 = plot(x,p,'k-', 'LineWidth', 1);
  xlabel('x (m)')
  ylabel('Pressure (kPa)')
  legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});

  set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10, 'XGrid', 'on', 'YGrid', 'on')
  
  % export figure with quality
  stmp = fullfile(figures_dir,'Mandel_pressure.png');
  exportgraphics(gcf,stmp,'Resolution',400)

  % Displacement DX
  figure('Position', [100, 100, 700, 700])
  hold on
  plotObj1 = plot(grid.coordinates(nodesX,1),dispXNum,'k.', 'LineWidth', 1, 'MarkerSize', 15);  
  plotObj2 = plot(x,ux,'k-', 'LineWidth', 1);
  xlabel('x (m)')
  ylabel('Displacement u_x (m)')
  ylim([0 7.e-5])
  legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'southeast');

  set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10, 'XGrid', 'on', 'YGrid', 'on')
  
  % export figure with quality
  stmp = fullfile(figures_dir,'Mandel_UX.png');
  exportgraphics(gcf,stmp,'Resolution',400)

  % Displacement DZ
  figure('Position', [100, 100, 700, 700])
  hold on
  plotObj1 = plot(dispZNum,grid.coordinates(nodesZ,3),'k.', 'LineWidth', 1, 'MarkerSize', 15);
  plotObj2 = plot(uz,z,'k-', 'LineWidth', 1);
  xlabel('Displacement u_z (m)')
  xlim([-10e-5 1e-5])
  ylabel('z (m)')
  legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});

  set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
  a = get(gca,'XTickLabel');
  set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10, 'XGrid', 'on', 'YGrid', 'on')

  % export figure with quality
  stmp = fullfile(figures_dir,'Mandel_UZ.png');
  exportgraphics(gcf,stmp,'Resolution',400)
end