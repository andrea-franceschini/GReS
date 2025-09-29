close all;
% clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

% isSPF = true;
isSPF = false;
% isXML = true;
isXML = false;

%% -------------------------- SET THE PHYSICS -------------------------
if isSPF
  physics = "SinglePhaseFlow";
  % model = ModelType("SinglePhaseFlow_FEM");
  model = ModelType("SinglePhaseFlow_FVTPFA");
else
  physics = "VariablySaturatedFlow";
  model = ModelType("VariabSatFlow_FVTPFA");
end

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = strcat(input_dir,'simParam.xml');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
fileName = strcat(input_dir,'Mesh/Column.msh');

% Import mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%% ----------------------------- MATERIALS -----------------------------
% Set the input file name
if isXML
  fileName = strcat(input_dir,'materials.xml');
else
  if isSPF
    fileName = strcat(input_dir,'SPF_materialsList.dat');
  else
    fileName = strcat(input_dir,'materialsList.dat');
  end
end

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%% ------------------------------ ELEMENTS -----------------------------
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
printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
  'folderName','Outputs');

%% ----------------------- Boundary Condition -----------------------------
% Creating and Appling boundaries conditions.
cond = struct('name',[],'type',[],'field',[],'values',[],'times',[]);
cond(1).name = 'Bottom';
cond(1).type = 'Dir';
cond(1).field = "latY0";%"bot";%"latY0";
% cond(1).times = 0.;
% cond(1).values = 1e4;%-10*9.8066e3;%1e4;%-10*9.8066e3;
cond(1).times = [0.;51840;77760.;129600.;259200.];
cond(1).values = [1e4;9e3;8e3;7e3;6e3;];%-10*9.8066e3;%1e4;%-10*9.8066e3;
cond(2).name = 'Top';
cond(2).type = 'Dir';
cond(2).field = "latYM";%"top";%"latYM";
cond(2).times = 0.;
cond(2).values = 1e3;%-0.75*9.8066e3;%1e3;%-0.75*9.8066e3;

fileName = setRichardsBC('Inputs',grid,cond,physics);
bound = Boundaries(fileName,model,grid);

%% ----------------------- Discretizer -----------------------------
% Create object handling construction of Jacobian and rhs of the model
domain = Discretizer('ModelType',model,...
  'SimulationParameters',simParam,...
  'DoFManager',dofmanager,...
  'Boundaries',bound,...
  'OutState',printUtils,...
  'Materials',mat,...
  'Grid',grid);

%% ----------------------- Initial Condition -----------------------------
% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = 1e3;%-10*9.8066e3;%1e4;%-10*9.8066e3;

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
postproc=false;
if postproc
  image_dir = fullfile(pwd,figures_dir);
  if ~isfolder(image_dir)
    mkdir(image_dir)
  end

  nrep = length(printUtils.results(:,:));
  nvars = length(printUtils.results(2,:).expPress);
  pressure = zeros(nvars,nrep);
  nvars = length(printUtils.results(2,:).expTime);
  t = zeros(nvars,nrep);
  for i=2:nrep
    pressure(:,i) = printUtils.results(i,:).expPress;
    t(:,i) = printUtils.results(i,:).expTime;
  end

  % Saving a temporary variabel.
  pressure = [printUtils.results.expPress];
  if ~isSPF
    saturation = [printUtils.results.expSat];
  end

  % Ajusting the time position.
  t = [printUtils.results.expTime];
  tind = 2:length(t);
  t_max = t(end);
  t = t(tind)/86400;
  tstr = strcat(num2str(t),' Days');

  %Getting pressure and saturation solution for specified time from MatFILE
  numbA = 5.;
  numbB = 5.;
  if isFEMBased(model,'Flow')
    tol = 1e-3;
    nodesP = find(abs(topology.coordinates(:,1)-numbA) < tol & abs(topology.coordinates(:,3)-numbB) < tol);
    [~,ind] = sort(topology.coordinates(nodesP,2));
    nodesP = nodesP(ind);

    pressplot = pressure(nodesP);

    % Values for normalized plots
    H = max(topology.coordinates(nodesP,2));

    % Location a column to be the plot position.
    pts = topology.coordinates(nodesP,2)/H;
  else
    tol = 0.4;
    nodesP = find(abs(topology.cellCentroid(:,1)-numbA) < tol & abs(topology.cellCentroid(:,3)-numbB) < tol);
    [~,ind] = sort(topology.cellCentroid(nodesP,2));
    nodesP = nodesP(ind);

    pressplot = pressure(nodesP);
    if ~isSPF
      swplot = saturation(nodesP);
    end

    % Values for normalized plots
    H = max(topology.cellCentroid(nodesP,2));

    % Location a column to be the plot position.
    pts = topology.cellCentroid(nodesP,2)/H;
  end
  weight = mat.getFluid().getFluidSpecWeight();





  %Plotting pressure head
  figure('Position', [100, 100, 700, 700])
  plot(pressplot/weight,ptsZ,'.-', 'LineWidth', 2, 'MarkerSize', 14);
  hold on
  xlabel('Head Pressure (m)')
  ylabel('Height (m)')
  legend(tstr, 'Location', 'northwest')
  set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
  % export figure with quality
  stmp = strcat(image_dir,'pressure','.png');
  exportgraphics(gcf,stmp,'Resolution',400)

  if ~isSPF
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
    stmp = strcat(image_dir, 'saturation', '.png');
    exportgraphics(gcf,stmp,'Resolution',400)
  end


  % To save the information necessary to compare the solution of GReS and MRST.
  % GReS_time = printUtils.results.expTime(tind,1);
  % GReS_pres = pressplot/weight;
  % GReS_satu = swplot;
  % GReS_elev = ptsZ;
  % save('Inputs/Solution/GReS_30.mat', 'GReS_time', 'GReS_pres', 'GReS_satu', 'GReS_elev');
end