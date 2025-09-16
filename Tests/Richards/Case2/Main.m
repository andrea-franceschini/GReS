close all;
% clear;
input_dir = 'Inputs';
output_dir = 'Outputs';
figures_dir = 'Figs';

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("VariabSatFlow_FVTPFA");

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = fullfile(input_dir,'simParam.dat');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
% fileName = fullfile(input_dir,'Mesh','Column.msh');
fileName = fullfile(input_dir,'Mesh','Column30.msh');

% Import mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%% ----------------------------- MATERIALS -----------------------------
% Set the input file name
fileName = fullfile(input_dir,'materialsList.dat');

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%% ------------------------------ ELEMENTS -----------------------------
% Define Gauss points
% GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
% elems = Elements(topology,GaussPts);
elems = Elements(topology);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model,topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

%% ----------------------- DOF Manager -----------------------------
% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

% Create and set the print utility
printUtils = OutState(model,topology,fullfile(input_dir,'outTime.dat'), ...
    'folderName','Outputs','flagMatFile',true);

%% ----------------------- Boundary Condition -----------------------------
% Creating and Appling boundaries conditions.
cond = struct('name',[],'type',[],'field',[],'values',[],'times',[]);
cond(1).name = 'Bottom';
cond(1).type = 'Dir';
cond(1).field = "bot";
cond(1).times = 0.;
cond(1).values = -10*9.8066e3;
cond(2).name = 'Top';
cond(2).type = 'Dir';
cond(2).field = "top";
cond(2).times = 0.;
cond(2).values = -0.75*9.8066e3;

fileName = setRichardsBC('Inputs',grid,cond);
bound = Boundaries(fileName,model,grid);

%% ----------------------- Discretizer -----------------------------
% Create object handling construction of Jacobian and rhs of the model
% linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
domain = Discretizer('ModelType',model,...
                     'SimulationParameters',simParam,...
                     'DoFManager',dofmanager,...
                     'Boundaries',bound,...
                     'OutState',printUtils,...
                     'Materials',mat,...
                     'Grid',grid);

%% ----------------------- Initial Condition -----------------------------
% Build a structure storing variable fields at each time step
% state = linSyst.setState();

% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = -10*9.8066e3;

%% ----------------------- Solver -----------------------------
% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(domain,'SaveRelError',true,'SaveBStepInf',true);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()

%% POST PROCESSING
postproc = true;
printFigs = true;
if postproc
  image_dir = fullfile(pwd,figures_dir);
    if ~isfolder(image_dir)
        mkdir(image_dir)
    end

    % Small modification - for the growning grid
    nrep = length(printUtils.results);
    nvars = length(printUtils.results(2).expPress);
    % nvars = length(printUtils.results(2).expSat);
    % nvars = length(printUtils.results(2).expTime);
    pressure = zeros(nvars,nrep);    
    saturation = zeros(nvars,nrep);    
    t = zeros(1,nrep);
    for i=2:nrep
       pressure(:,i) = printUtils.results(i).expPress;
       saturation(:,i) = printUtils.results(i).expSat;
       t(i) = printUtils.results(i).expTime;
    end

    % Ajusting the time position.
    % t = printUtils.results.expTime;
    tind = 2:length(t);
    t_max = t(end);
    t = t(tind)/86400;
    tstr = strcat(num2str(t'),' Days');

    %Getting pressure and saturation solution for specified time from MatFILE
    numb = 0.;
    tol = 0.01;
    % nodesP = find(abs(elems.cellCentroid(:,1)-numb) < tol & abs(elems.cellCentroid(:,2)-numb) < tol);
    nodesP = find(abs(elems.mesh.cellCentroid(:,1)-numb) < tol & abs(elems.mesh.cellCentroid(:,2)-numb) < tol);
    pressplot = pressure(nodesP,tind);
    swplot = saturation(nodesP,tind);

    % Values for normalized plots
    H = max(topology.coordinates(:,3));
    weight = mat.getFluid().getFluidSpecWeight();

    % Location a column to be the plot position.
    ptsZ = elems.mesh.cellCentroid(nodesP,3);

    if printFigs
      %Plotting pressure head
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

      %Plotting saturation
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

    % To save the information necessary to compare the solution of GReS and MRST.
    % GReS_time = [printUtils.results(tind).expTime];
    % GReS_pres = pressplot/weight;
    % GReS_satu = swplot;
    % GReS_elev = ptsZ;
    % save('Inputs/Solution/GReS_30.mat', 'GReS_time', 'GReS_pres', 'GReS_satu', 'GReS_elev');
end

% run("Validation.m")