close all;
clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("Poromechanics_FEM");

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = strcat(input_dir,'simParam.dat');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
fileName = strcat(input_dir,'Mesh/','Beam.msh');
% fileName = strcat(input_dir,'Mesh/','Beam_20_3.msh');

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
% cond = struct('name',[],'type',[],'physics',[],'field',[],'values',[],'times',[]);
% cond(1).name = 'prescribe';
% cond(1).type = 'Dir';
% cond(1).physics = 'Poro';
% cond(1).field = "latX0";
% cond(1).times = 0.;
% cond(1).values = 0.;
% 
% cond(2).name = 'desloc';
% cond(2).type = 'Dir';
% cond(2).physics = 'Poro';
% cond(2).field = "latXM";
% cond(2).times = 0.;
% cond(2).values = -10.;
% 
% fileName = setBC_Nodal(input_dir,grid,cond);

fileName = ["Inputs/Boundary_beam_10/force.dat","Inputs/Boundary_beam_10/prescribe.dat"];
% fileName = ["Inputs/Boundary_beam_20_3/force.dat","Inputs/Boundary_beam_20_3/prescribe.dat"];
% fileName = ["Inputs/Boundary_beam_20_3/desloc.dat","Inputs/Boundary_beam_20_3/prescribe.dat"];
% fileName = ["Inputs/Boundary_beam_20_3/desloc_node.dat","Inputs/Boundary_beam_20_3/prescribe.dat"];

bound = Boundaries(fileName,model,grid);

%% -------------------------- PREPROCESSING ----------------------------
% Create and set the print utility
printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
    'folderName','Outputs');

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
if true
    expDispl = printUtils.results.expDispl;
    expTime = printUtils.results.expTime;
    elems = [51;52];
    nodes = [1;9;10;11;12;13;14;15;16;17;2];
    % elems = [430;431;432;433;434;435;436;437;438];
    % nodes = [1;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;2];
    % [~, order,~] = unique(nodes);
    dispPlot = expDispl(nodes*3,:);

    figure()
    plot(dispPlot(:,2),'-ko','LineWidth',1,'MarkerSize', 5)
    xlabel('Time (days)')
    ylabel('Vertical displacements (mm)')
    xlim([0 1.05])
    ylim([-40 5])
    set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
    % export figure
    stmp = strcat(figures_dir, 'SurfLoad_dispTime', '.png');
    exportgraphics(gcf,stmp,'Resolution',400)
end