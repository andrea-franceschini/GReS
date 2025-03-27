close all;
clear;
input_dir = 'Inputs/';
output_dir = 'Outputs/';

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("VariabSatFlow_FVTPFA");

%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = strcat(input_dir,'simParam.dat');
simParam = SimulationParameters(fileName,model);

%% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();

% Set the input file name
% fileName = strcat(input_dir,'Column.msh');
fileName = strcat(input_dir,'Mesh/Column160.msh');

% Import mesh data into the Mesh object
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
state.pressure(:) = -10*9.8066e3;

% Create and set the print utility
printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
    'folderName','Outputs');

printState(printUtils,state)

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

Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()

%% POST PROCESSING
postproc=true;
if postproc
    image_dir = strcat(pwd,'/',output_dir,'Images');
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

    % Location a column to be the plot position.

    %Getting pressure and saturation solution for specified time from MatFILE
    numb = 0.;
    tol = 0.01;
    nodesP = find(abs(elems.cellCentroid(:,1)-numb) < tol & abs(elems.cellCentroid(:,2)-numb) < tol);
    pressplot = pressure(nodesP,tind);
    swplot = saturation(nodesP,tind);

    % Values for normalized plots
    H = max(topology.coordinates(:,3));

    %Plotting solution
    ptsZ = topology.coordinates(nodesP,3);

    figure(1)
    plot(-pressplot,ptsZ/H,'.-', 'LineWidth', 2, 'MarkerSize', 14);
    hold on
    xlabel('p/p_{top}')
    ylabel('z/H')
    legend(tstr)
    % grid on
    set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 16);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 14)
    % export figure with quality
    stmp = strcat(output_dir,'Images/','Varelha_head_pressure','.png');
    exportgraphics(gcf,stmp,'Resolution',400)

end

if false
    %Post processing using MAT-FILE

    %Plotting solution
    if isFVTPFABased(model,'Flow')
        ptsY = elems.cellCentroid(nodesP,3);
    else
        ptsY = topology.coordinates(nodesP,3);
    end
    figure(1)
    plot(-pressplot,ptsY/H,'.-', 'LineWidth', 1, 'MarkerSize', 10);
    hold on
    xlabel('p/p_{top}')
    ylabel('z/H')
    legend(tstr)
    % grid on
    set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 12)
    % export figure with quality
    stmp = strcat(output_dir,'Images/','Richards_pressure_old','.png');
    exportgraphics(gcf,stmp,'Resolution',400)

    figure(2)
    plot(swplot,ptsY/H,'.-', 'LineWidth', 1, 'MarkerSize', 10);
    hold on
    xlabel('Saturation S_w')
    ylabel('z/H')
    legend(tstr, 'Location', 'southwest')
    str = strcat('t = ',tstr);
    % grid on
    set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 12)
    % export figure with quality
    stmp = strcat(output_dir,'Images/', 'Richards_staturation_old', '.png');
    exportgraphics(gcf,stmp,'Resolution',400)

end