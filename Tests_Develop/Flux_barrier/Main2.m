close all;
clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("SinglePhaseFlow_FVTPFA");
% model = ModelType("SinglePhaseFlow_FEM");

%% ----------------------- SIMULATION PARAMETERS ----------------------
 fileName = strcat(input_dir,'simParam.dat');
% fileName = strcat(input_dir,'simParam.xml');
% fileName = strcat(input_dir,'simParam2.xml');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
%fileName = strcat(input_dir,'Mesh/Fault.msh');
%fileName = strcat(input_dir,'Mesh/Pillar.msh');
fileName = strcat(input_dir,'Mesh/Double_Pillar.msh');
%fileName = strcat(input_dir,'Mesh/Double_Pillar_2.msh');

% Import mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%% ----------------------------- MATERIALS -----------------------------
% Set the input file name
fileName = strcat(input_dir,'materialsList.dat');
% fileName = strcat(input_dir,'materials.xml');

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
physics = "SinglePhaseFlow";
% physics = "VariablySaturatedFlow";
cond = struct('name',[],'type',[],'field',[],'values',[],'times',[]);
cond(1).name = 'BoundA';
cond(1).type = 'Dir';
cond(1).field = "X0";
cond(1).times = 0.;
cond(1).values = 10.;
cond(2).name = 'BoundB';
cond(2).type = 'Dir';
cond(2).field = "Xm";
cond(2).times = 0.;
cond(2).values = 0.;

%fileName = ["BC/Time_BoundA_0.dat","BC/Time_BoundB_0.dat"];
fileName = setBoundaryC('Inputs',grid,cond,physics);
bound = Boundaries(fileName,model,grid);

%% ----------------------- Discretizer -----------------------------
% Create object handling construction of Jacobian and rhs of the model

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
% set initial conditions directly modifying the state object
domain.state.data.pressure(:) = 0;

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

      % Values for normalized plots
      H = max(topology.cellCentroid(nodesP,2));

      % Location a column to be the plot position.
      pts = topology.cellCentroid(nodesP,2)/H;
    end
    weight = mat.getFluid().getFluidSpecWeight();

    %Plotting pressure head
    figure(1)
    plot(-pressplot/weight,pts,'.-', 'LineWidth', 2, 'MarkerSize', 14);
    hold on
    xlabel('pressure (m)')
    ylabel('height (m)')
    legend(tstr)
    set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
    % export figure with quality
    stmp = strcat(image_dir,'Varelha_head_pressure','.png');
    exportgraphics(gcf,stmp,'Resolution',400)
end
