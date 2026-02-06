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
topology = Mesh();

% Import mesh data into the Mesh object
topology.importMesh(fullfile(input_dir,'Mesh',"Column30.msh"));

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
printUtils = OutState(fullfile(input_dir,'output.xml'));

% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('Grid',grid,...
                     'Materials',mat,...
                     'Boundaries',bound);

domain.addPhysicsSolver(fullfile(input_dir,'solver.xml'));

% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = -9.8066e4;

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
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

    nrep = length(printUtils.matFile);
    nvars = length(printUtils.matFile(2).pressure);
    pressure = zeros(nvars,nrep);    
    saturation = zeros(nvars,nrep);    
    t = zeros(1,nrep);
    for i=1:nrep
       pressure(:,i) = printUtils.matFile(i).pressure;
       saturation(:,i) = printUtils.matFile(i).saturation;
       t(i) = printUtils.matFile(i).time;
    end

    % Adjusting the time position.
    tind = 1:length(t);
    t_max = t(end);
    t = t(tind)/86400;
    tstr = strcat(num2str(t'),' Days');

    % Getting pressure and saturation solution for specified time from MatFILE
    numb = 0.;
    tol = 0.01;
    nodesP = find(abs(elems.mesh.cellCentroid(:,1)-numb) < tol & abs(elems.mesh.cellCentroid(:,2)-numb) < tol);
    pressplot = pressure(nodesP,tind);
    swplot = saturation(nodesP,tind);

    % Values for normalized plots
    H = max(topology.coordinates(:,3));
    weight = mat.getFluid().getSpecificWeight();

    % Location a column to be the plot position.
    ptsZ = elems.mesh.cellCentroid(nodesP,3);

    if printFigs
      % Plotting pressure head
      figure('Position', [100, 100, 700, 700])
      plot(pressplot/weight,ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
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
      plot(swplot,ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
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