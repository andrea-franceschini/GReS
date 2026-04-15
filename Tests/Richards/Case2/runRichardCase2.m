close all;
% clear;
output_dir = 'Output';
input_dir = 'Input';
figures_dir = fullfile(output_dir,"Images");

%% ------------------------------ Set up the Domain -----------------------
% Set the simulation parameters for the non-linear solver.
simParam = SimulationParameters(fullfile(input_dir,'simparam.xml'));

% Create an object of the Materials class and read the materials file
mat = Materials(fullfile(input_dir,'materials.xml'));

% Create the Mesh object
grid = Grid();

% Import mesh data into the Mesh object
grid.importMesh(fullfile(input_dir,'Mesh',"Column30.msh"));


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

% Set initial conditions directly modifying the state object
domain.state.data.pressure(:) = -9.8066e4;

% Set and solve the simulation
solver = NonLinearImplicit('simulationparameters',simParam,...
                           'domains',domain,...
                           'output',printUtils);
solver.simulationLoop();

%% --------------------- Post Processing the Results ----------------------
postproc = true;
printFigs = true;
if postproc
  image_dir = fullfile(pwd,figures_dir);
    if ~isfolder(image_dir)
        mkdir(image_dir)
    end

    % Getting pressure and saturation solution for specified time from MatFILE
    numb = 0.;
    tol = 0.01;
    nodesP = find(abs(elems.mesh.cellCentroid(:,1)-numb) < tol & abs(elems.mesh.cellCentroid(:,2)-numb) < tol);

    tstr = strcat(num2str((printUtils.timeList/86400)'),' Days');
    nrep = length(printUtils.timeList);
    nvars = length(nodesP);
    press = zeros(nvars,nrep);
    sw = zeros(nvars,nrep);
    for i=1:length(printUtils.timeList)
      press(:,i) = printUtils.results(i).pressure(nodesP);
      sw(:,i) = printUtils.results(i).saturation(nodesP);
    end

    % Values for normalized plots
    % H = max(topology.coordinates(:,3));
    weight = mat.getFluid().getSpecificWeight();

    % Location a column to be the plot position.
    ptsZ = elems.mesh.cellCentroid(nodesP,3);

    if printFigs
      % Plotting pressure head
      figure('Position', [100, 100, 700, 700])
      plot(press/weight,ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
      hold on
      xlabel('Head Pressure (m)')
      ylabel('Height (m)')
      legend(tstr, 'Location', 'northwest')
      set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
      % export figure with quality
      stmp = fullfile(image_dir,'pressure.png');
      exportgraphics(gcf,stmp,'Resolution',400)

      % Plotting saturation
      figure('Position', [100, 100, 700, 700])
      plot(sw,ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
      hold on
      xlabel('Saturation')
      ylabel('Height (m)')
      legend(tstr, 'Location', 'northwest')
      str = strcat('t = ',tstr);
      set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
      % export figure with quality
      stmp = fullfile(image_dir,'saturation.png');
      exportgraphics(gcf,stmp,'Resolution',400)
    end
end