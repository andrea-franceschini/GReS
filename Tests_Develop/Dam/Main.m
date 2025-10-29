close all;
% clear;
input_dir = 'Inputs';
output_dir = 'Outputs';
figures_dir = 'Figs';

isSPF = true;

%% -------------------------- SET THE PHYSICS -------------------------
if isSPF
  model = ModelType("SinglePhaseFlow_FVTPFA");
else
  model = ModelType("VariabSatFlow_FVTPFA");
end

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = fullfile(input_dir,'simParam.xml');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
fileName = fullfile(input_dir,'Mesh','dam.msh');

% Import mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%% ----------------------------- MATERIALS -----------------------------
% Set the input file name
if isSPF
  fileName = fullfile(input_dir,'materialsListSPF.dat');
else
  fileName = fullfile(input_dir,'materialsListVSF.dat');
end

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%% ------------------------------ ELEMENTS -----------------------------
% Define Gauss points
% GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

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

% cond(1).name = 'wLev';
% cond(1).type = 'Dir';
% cond(1).field = "topA";
% cond(1).times = 0.;
% % cond(1).values = -0.75*9.8066e3;
% cond(1).values = 1e3;

% cond(1).name = 'wLev';
% cond(1).type = 'Spg';
% cond(1).field = "topA";
% cond(1).times = 0.;
% % cond(1).values = -0.75*9.8066e3;
% cond(1).values = 23.5;

% cond(2).name = 'Ux0';
% cond(2).type = 'Dir';
% cond(2).field = "latX0";
% cond(2).times = 0.;
% cond(2).values = -10*9.8066e3;
% cond(2).name = 'Ux0';
% cond(2).type = 'Spg';
% cond(2).field = "latX0";
% cond(2).times = 0.;
% cond(2).values = 20.;

% cond(3).name = 'Uz0';
% cond(3).type = 'Neu';
% cond(3).field = "bot";
% cond(3).times = 0.;
% cond(3).values = 0;
% 
% cond(3).name = 'NeuB';
% cond(3).type = 'Neu';
% cond(3).field = "topB";
% cond(3).times = 0.;
% cond(3).values = 0.;
% 
% cond(4).name = 'NeuC';
% cond(4).type = 'Neu';
% cond(4).field = "topC";
% cond(4).times = 0.;
% cond(4).values = 0.;


% cond(1).name = 'Ux0';
% cond(1).type = 'Dir';
% cond(1).field = "latX0";
% cond(1).times = 0.;
% cond(1).values = 1e6;

% cond(2).name = 'UxM';
% cond(2).type = 'Dir';
% cond(2).field = "latXM";
% cond(2).times = 0.;
% cond(2).values = 1e5;


% cond(1).name = 'Ux0';
% cond(1).type = 'Spg';
% cond(1).field = "latX0";
% cond(1).times = 0.;
% cond(1).values = 20.;

cond(1).name = 'wLev';
cond(1).type = 'Spg';
cond(1).field = "topA";
cond(1).times = 0.;
cond(1).values = 4.5;


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
% domain.state.data.pressure(:) = -10*9.8066e3;
% domain.state.data.pressure(:) = 1e5;
domain.state.data.pressure(:) = 0.;

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