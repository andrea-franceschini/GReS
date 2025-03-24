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
fileName = strcat(input_dir,'Mesh/Column.msh');

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
z = elems.cellCentroid(:,3);
gamma_w = getFluid(mat).getFluidSpecWeight();
wLev = 9.; % level of the water table
state.pressure = gamma_w*(wLev-z);

% Create and set the print utility
printUtils = OutState(model,topology,strcat(input_dir,'outTime.dat'), ...
    'folderName','Outputs');

printState(printUtils,state)

% Creating and Appling boundaries conditions.
cond = struct('name',[],'type',[],'field',[],'values',[],'times',[]);
% % Old mesh file.
% cond(1).name = 'Bottom';
% cond(1).type = 'Dir';
% cond(1).field = "lFace";
% cond(1).times = [0. 5. 10.];
% cond(1).values = [87. 0. 0.];
% cond(2).name = 'Top';
% cond(2).type = 'Dir';
% cond(2).field = "lFace";
% cond(2).times = [0. 5. 10.];
% cond(2).values = [87. 0. 0.];

% % New gmsh file.
cond(1).name = 'Bottom';
cond(1).type = 'Dir';
cond(1).field = "bot";
cond(1).times = [0. 5. 10.];
cond(1).values = [87. 0. 0.];
% cond(2).name = 'Countors';
% cond(2).type = 'Neu';
% cond(2).field = ["latXM","latYM","latX0","latY0"];
% cond(2).times = 0.;
% cond(2).values = 0.;
% cond(3).name = 'Top';
% cond(3).type = 'Dir';
% cond(3).field = "top";
% cond(3).times = 0.;
% cond(3).values = 0.;
% cond(4).name = 'Spg';
% cond(4).type = 'Spg';
% cond(4).field = ["latXM","latYM","latX0","latY0"];
% cond(4).times = [0. 5. 10.];
% cond(4).values = [87. 0. 0.];

fileName = setRichardsBC('Inputs',grid,cond);
bound = Boundaries(fileName,model,grid);

Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);

% Solve the problem
[simState] = Solver.NonLinearLoop();

% Finalize the print utility
printUtils.finalize()

%% POST PROCESSING

image_dir = strcat(pwd,'/',output_dir,'Images');
if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir(image_dir)
else
    mkdir(image_dir)
end

%Post processing using MAT-FILE 

% elem vector containing elements centroid along vertical axis
if isFEMBased(model,'Flow')
    nodesP = nodesU;
else
    numb = 0.125;
    tol = 0.01;
    nodesP = find(abs(elems.cellCentroid(:,1)-numb) < tol & abs(elems.cellCentroid(:,2)-numb) < tol);
    [~,ind] = sort(elems.cellCentroid(nodesP,3));
    nodesP = nodesP(ind);
end

press = printUtils.results.expPress;
sw = printUtils.results.expSat;
t = printUtils.results.expTime;
tind = 2:length(t);
t_max = t(end);
t = t(tind)/t_max;


tstr = strcat(num2str(t),' T');
%Getting pressure and saturation solution for specified time from MatFILE
pressplot = press(nodesP,tind);
swplot = sw(nodesP,tind);

% Values for normalized plots
H = max(topology.coordinates(:,3));

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