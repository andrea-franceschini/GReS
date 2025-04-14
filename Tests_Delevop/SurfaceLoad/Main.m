close all;
clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = strcat(input_dir,'simParam.dat');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
% fileName = strcat(input_dir,'Mesh/','Test01_hexa.msh');
fileName = strcat(input_dir,'Mesh/','QrCircle.msh');

% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

%% ----------------------------- MATERIALS -----------------------------
% Set the input file name
fileName = strcat(input_dir,'materialsList.dat');

% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

%% ------------------------------ ELEMENTS -----------------------------
% Define Gauss points
GaussPts = Gauss(12,2,3);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);

% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);

% Degree of freedom manager 
dofmanager = DoFManager(topology,model);
% dofmanager = DoFManager_new(topology,model,'dof.dat');

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

% Build a structure storing variable fields at each time step
resState = linSyst.setState();

%% ------------------------ BOUNDARY CONDITIONS ------------------------
% Set the input file
% Creating and Appling boundaries conditions.
cond = struct('name',[],'type',[],'physics',[],'field',[],'values',[],'times',[]);
cond(1).name = 'Bot';
cond(1).type = 'Dir';
cond(1).physics = 'Poro';
cond(1).field = "bot";
cond(1).times = 0.;
cond(1).values = 0.;

cond(2).name = 'TopI';
cond(2).type = 'Dir';
cond(2).physics = 'Poro';
cond(2).field = "topI";
cond(2).times = [0. 3. 100.];
cond(2).values = [0. -10. -10.];

cond(3).name = 'TopE';
cond(3).type = 'Dir';
cond(3).physics = 'Poro';
cond(3).field = "topE";
cond(3).times = 0.;
cond(3).values = 0.;

cond(4).name = 'LateralX';
cond(4).type = 'Neu';
cond(4).physics = 'Flow';
cond(4).field = "latX";
cond(4).times = 0.;
cond(4).values = 0.;

cond(5).name = 'LateralY';
cond(5).type = 'Neu';
cond(5).physics = 'Flow';
cond(5).field = "latY";
cond(5).times = 0.;
cond(5).values = 0.;

cond(6).name = 'LateralC';
cond(6).type = 'Neu';
cond(6).physics = 'Flow';
cond(6).field = "latC";
cond(6).times = 0.;
cond(6).values = 0.;

fileName = setBC(input_dir,grid,cond);

% fileName = ["dir_Poro_fixed.dat","neuSurf_Poro_TopLoad.dat","dir_Poro_Sym.dat","dir_Flow_TopDrained.dat"];
%
bound = Boundaries(fileName,model,grid);
% bound = Boundaries(fileName,model,grid,dofmanager);

%% -------------------------- PREPROCESSING ----------------------------
% Create and set the print utility
printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
    'folderName','Outputs');
% printUtils = OutState(model,mat,grid,'outTime.dat');

% Print the reservoir initial state
printUtils.printState(resState);

%% ---------------------------- SOLUTION -------------------------------
% Create the object handling the (nonlinear) solution of the problem
NSolv = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,resState,linSyst,'GaussPts',GaussPts,'SaveRelError',true,'SaveBStepInf',true);

% Solve the problem
[simState] = NSolv.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()


%% -------------------------- POST PROCESSING ------------------------------
if false
    expPress = printUtils.results.expPress;
    expDispl = printUtils.results.expDispl;
    expTime = printUtils.results.expTime;
    elems = [3427; 2467; 1507; 547];
    nodes = [359; 287; 215; 143];
    pressplot = expPress(elems,:);
    dispPlot = expDispl(nodes*3,:);
    figure(1)
    plot(expTime,pressplot,'-ko','LineWidth',1,'MarkerSize', 5)
    xlabel('Time (days)')
    ylabel('Pressure (kPa)')
    xlim([0 10.5])
    ylim([-1 6])
    grid on
    set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
    % export figure with quality
    stmp = strcat(figures_dir, 'SurfLoad_pressure', '.png');
    exportgraphics(gcf,stmp,'Resolution',400)

    figure(2)
    plot(expTime,1000*dispPlot,'-ko','LineWidth',1,'MarkerSize', 5)
    xlabel('Time (days)')
    ylabel('Vertical displacements (mm)')
    xlim([0 10.5])
    ylim([-40 5])
    grid on
    set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
    % export figure
    stmp = strcat(figures_dir, 'SurfLoad_dispTime', '.png');
    exportgraphics(gcf,stmp,'Resolution',400)



    figure(3)
    timesInd = [3 5 7 17];
    nodes_subsidence = intersect(find(topology.coordinates(:,3) ==50),find(topology.coordinates(:,2) ==0));
    x = sort(topology.coordinates(nodes_subsidence,1));
    subs = sort(expDispl(3*nodes_subsidence,timesInd),'ascend');
    set(0,'DefaultAxesColorOrder',[0 0 0],...
        'DefaultAxesLineStyleOrder','-|--|:|-.')
    plot(x,1000*subs, 'LineWidth', 1);
    xlabel('r (m)')
    ylabel('Vertical displacements (mm)')
    xlim([0 50])
    ylim([-60 5])
    time_string = "Day  " + expTime(timesInd);
    legend(time_string, 'Location', 'southeast')
    grid on
    set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
    % export figure with quality
    stmp = strcat(figures_dir, 'SurfLoad_dispR', '.png');
    exportgraphics(gcf,stmp,'Resolution',400)

end