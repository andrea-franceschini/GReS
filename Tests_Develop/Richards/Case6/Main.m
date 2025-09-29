close all;
% clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

% isSPF = true;
isSPF = false;

%% -------------------------- SET THE PHYSICS -------------------------
if isSPF
   model = ModelType("SinglePhaseFlow_FVTPFA");
   % model = ModelType("SinglePhaseFlow_FEM");
else
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
if isSPF
   fileName = strcat(input_dir,'SPF_materialsList.dat');
else
   fileName = strcat(input_dir,'materialsList.dat');
end

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%% ------------------------------ ELEMENTS -----------------------------
% Define Gauss points
GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model,topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

% Build a structure storing variable fields at each time step
state = linSyst.setState();

% set initial conditions directly modifying the state object
state.pressure(:) = 1e3;%-10*9.8066e3;%1e4;%-10*9.8066e3;

% Create and set the print utility
printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
    'folderName','Outputs');

printState(printUtils,state)

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

fileName = setRichardsBC('Inputs',grid,cond);
bound = Boundaries(fileName,model,grid);

Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,...
   printUtils,state,linSyst,'GaussPts',GaussPts,'SaveConverg',true);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()

%% POST PROCESSING
postproc=false;
if postproc
    image_dir = strcat(pwd,'/',figures_dir);
    if ~isfolder(image_dir)
        mkdir(image_dir)
    end

    % Saving a temporary variabel.
    pressure = printUtils.results.expPress;
    saturation = printUtils.results.expSat;

    % Ajusting the time position.
    t = printUtils.results.expTime;
    tind = 2:length(t);
    t_max = t(end);
    t = t(tind)/86400;
    tstr = strcat(num2str(t),' Days');

    %Getting pressure and saturation solution for specified time from MatFILE
    numb = 0.;
    tol = 0.01;
    nodesP = find(abs(elems.cellCentroid(:,1)-numb) < tol & abs(elems.cellCentroid(:,2)-numb) < tol);
    pressplot = pressure(nodesP,tind);
    swplot = saturation(nodesP,tind);

    % Values for normalized plots
    H = max(topology.coordinates(:,3));
    weight = mat.getFluid().getFluidSpecWeight();

    % Location a column to be the plot position.
    ptsZ = elems.cellCentroid(nodesP,3);

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

    % To save the information necessary to compare the solution of GReS and MRST.
    % GReS_time = printUtils.results.expTime(tind,1);
    % GReS_pres = pressplot/weight;
    % GReS_satu = swplot;
    % GReS_elev = ptsZ;
    % save('Inputs/Solution/GReS_30.mat', 'GReS_time', 'GReS_pres', 'GReS_satu', 'GReS_elev');
end