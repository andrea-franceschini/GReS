close all;
clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';
figures_dir = 'Figs/';

%% -------------------------- SET THE PHYSICS -----------------------------
model = ModelType("VariabSatFlow_FVTPFA");

%% ----------------------- SIMULATION PARAMETERS --------------------------
fileName = strcat(input_dir,'simParamMesh.dat');
simParam = SimulationParameters(fileName,model);

%% ----------------------------- MATERIALS --------------------------------
% Set the input file name
fileMaterial = strcat(input_dir,'materialsList.dat');

%% ------------------------------ ELEMENTS --------------------------------
% Define Gauss points
GaussPts = Gauss(12,2,3);

%% ------------------------------ BOUNDARY CONDITIONS ---------------------
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

%% ------------------------------ MESH ------------------------------------
meshDir = 'Mesh/';
meshName = 'Column';
meshList = [10 20 40 80 160];

% Structure to save the results for the convergence
result = struct('pressure',[],'saturation',[],'height',[]);

for i=1:length(meshList)
    % Create the Mesh object
    topology = Mesh();

    % Set the input file name
    fileName = strcat(input_dir,meshDir,meshName,int2str(meshList(i)),'.msh');

    % Import mesh data into the Mesh object
    topology.importGMSHmesh(fileName);

    % Create an object of the Materials class and read the materials file
    mat = Materials(model,fileMaterial);

    % Create an object of the "Elements" class and process the element properties
    elems = Elements(topology,GaussPts);

    % Create an object of the "Faces" class and process the face properties
    faces = Faces(model,topology);

    % Wrap Mesh, Elements and Faces objects in a structure
    grid = struct('topology',topology,'cells',elems,'faces',faces);

    % Degree of freedom manager
    dofmanager = DoFManager(topology,model);

    % Create object handling construction of Jacobian and rhs of the model
    linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

    % Build a structure storing variable fields at each time step
    state = linSyst.setState();

    % set initial conditions directly modifying the state object
    state.pressure(:) = -10*9.8066e3;

    % Create and set the print utility
    printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
        'folderName','Outputs','writeVtk',false);

    printState(printUtils,state)

    % Appling boundaries conditions.
    fileName = setRichardsBC('Inputs',grid,cond);
    bound = Boundaries(fileName,model,grid);

    % Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);
    Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,'GaussPts',GaussPts,'SaveRelError',true);

    % Solve the problem
    [simState] = Solver.NonLinearLoop();

    % Finalize the print utility
    printUtils.finalize()

    % Retriving the analisys information
    pressure = printUtils.results.expPress;
    saturation = printUtils.results.expSat;

    %Getting pressure and saturation solution for specified time from MatFILE
    numb = 0.;
    tol = 0.01;
    nodesP = find(abs(elems.cellCentroid(:,1)-numb) < tol & abs(elems.cellCentroid(:,2)-numb) < tol);
    result(i).pressure = pressure(nodesP,2);
    result(i).saturation = saturation(nodesP,2);
    result(i).height = elems.cellCentroid(nodesP,3);
end

%% POST PROCESSING
image_dir = strcat(pwd,'/',figures_dir);
if ~isfolder(image_dir)
    mkdir(image_dir)
end

% Ajusting the legend
tstr = strcat(int2str(meshList(:)),' Cells');

% Ajusting the pressure
weight = mat.getFluid().getFluidSpecWeight();

%Plotting pressure head
figure('Position', [100, 100, 700, 700])
hold on
for i=1:length(result)
    plot(result(i).pressure/weight,result(i).height,'.-', 'LineWidth', 2, 'MarkerSize', 14);
end
xlabel('Pressure (m)')
ylabel('Height (m)')
legend(tstr, 'Location', 'northwest')
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(figures_dir,'mesh_pressure','.png');
exportgraphics(gcf,stmp,'Resolution',400)

%Plotting saturation
figure('Position', [100, 100, 700, 700])
hold on
for i=1:length(result)
    plot(result(i).saturation,result(i).height,'.-', 'LineWidth', 2, 'MarkerSize', 14);
end
xlabel('Saturation')
ylabel('Height (m)')
legend(tstr, 'Location', 'northwest')
str = strcat('t = ',tstr);
set(gca,'FontName', 'Liberation Serif', 'FontSize', 16, 'XGrid', 'on', 'YGrid', 'on')
% export figure with quality
stmp = strcat(figures_dir, 'mesh_saturation', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

% To save the information necessary to compare the solution of GReS and MRST.
% GReS_mesh = meshList;
% GReS_result = result;
% GReS_weight = weight;
% save('Inputs/Solution/GReS_Mesh_Conv.mat', 'GReS_mesh', 'GReS_result', 'GReS_weight');
