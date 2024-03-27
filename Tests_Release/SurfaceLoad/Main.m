close all;
clear;

warning('off','MATLAB:nearlySingularMatrix');

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "simParam.dat";
simParam = SimulationParameters(model,fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'Test01_hexa.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'materialsList.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);
%
%------------------------------ ELEMENTS -----------------------------
%
% Define Gauss points
GaussPts = Gauss(12,2,3);
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
dofmanager = DoFManager(topology,model);
%
%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["dir_Poro_fixed.dat","neuSurf_Poro_TopLoad.dat","dir_Poro_Sym.dat","dir_Flow_TopDrained.dat"];
%
bound = Boundaries(fileName,model,grid,dofmanager);
%
%-------------------------- PREPROCESSING ----------------------------
%
resState = State(model,grid,mat,GaussPts);
% Create and set the print utility
printUtils = OutState(model,mat,grid,'outTime.dat');
%
% Print the reservoir initial state
printUtils.printState(resState);
%
% ---------------------------- SOLUTION -------------------------------
%

% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,resState,GaussPts);
%
% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%

%%
% -------------------------- POST PROCESSING ------------------------------

image_dir = strcat(pwd,'/Images');
if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir Images
else
    mkdir Images
end

expPress = printUtils.m.expPress;
expDispl = printUtils.m.expDispl;
expTime = printUtils.m.expTime;
elems = [3427; 2467; 1507; 547];
nodes = [359; 287; 215; 143];
pressplot = expPress(elems,:);
dispPlot = expDispl(nodes*3,:);
%
%
figure(1)
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
colororder(newcolors)
plot(expTime,pressplot,'-o','LineWidth',1,'MarkerSize', 5)
xlabel('Time (days)')
ylabel('Pressure (kPa)')
xlim([0 10.5])
ylim([-1 6])
legend('Upper Silty Clay', 'Sand', 'Lower Silty Clay', 'Silt')
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'SurfLoad_pressure', '.png');
exportgraphics(gcf,stmp,'Resolution',400)
%
%
figure(2)
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
colororder(newcolors)
plot(expTime,1000*dispPlot,'-o','LineWidth',1,'MarkerSize', 5)
xlabel('Time (days)')
ylabel('Vertical displacements (mm)')
xlim([0 10.5])
ylim([-40 5])
legend('Upper Silty Clay', 'Sand', 'Lower Silty Clay', 'Silt')
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure
stmp = strcat('Images\', 'SurfLoad_dispTime', '.png');
exportgraphics(gcf,stmp,'Resolution',400)
%
%
figure(3)
newcolors = 'k';
colororder(newcolors)
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
stmp = strcat('Images\', 'SurfLoad_dispR', '.png');
exportgraphics(gcf,stmp,'Resolution',400)


% figure(4)
% t = [0 1 6 7 100];
% q = [0 -90000 -90000 0 0];
% plot(t,q,'.-','LineWidth',1,'MarkerSize',15)
% xlim([-1 10])
% ylim([-100000 5000])
% set(gca, 'YDir','reverse')
% xlabel('Pressione (kPa)')
% ylabel('Flusso prescritto (m^3/anno)')

% figure(5)
% t = [0 3 100];
% q = [0 10 10];
% plot(t,q,'.-','LineWidth',1,'MarkerSize',15)
% xlim([-1 10])
% ylim([-1 12])
% xlabel('Tempo (giorni)')
% ylabel('Carico agente (kPa)')
