close all;
clear;

warning('off','MATLAB:nearlySingularMatrix');

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType(["SinglePhaseFlow_FEM","Poromechanics_FEM"]);
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "simParam.dat";
simParam = SimulationParameters(fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'ReservoirTest_Hexa.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'materialsListElastic.dat';
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
%----------------------------- DOF MANAGER -----------------------------
fileName = 'dof.dat';
if strcmp(fileName,'dof.dat')
    dofmanager = DoFManager(topology, model, fileName);
else
    dofmanager = DoFManager(topology, model);
end

%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
if strcmp(fileName, 'dof.dat')
    fileName = ["bottom_fixed.dat","flux.dat","lateral_fix.dat","Impermeable.dat"];
else
    fileName = ["bottom_fixed.dat","flux.dat","lateral_fix.dat","Impermeable.dat", "top_drained.dat"];
end
%
bound = Boundaries(fileName,model,grid,dofmanager);

% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat,GaussPts);

% Create and set the print utility
printUtils = OutState(model,mat,grid,'outTime.dat','printOn','out');
%
% Print the reservoir initial state
printUtils.printState(resState);
%
% ---------------------------- SOLUTION -------------------------------
%
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,resState,linSyst,GaussPts);
%
% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()

%% POST PROCESSING

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

%find nodes in vertical symmetry axis
tmp1=topology.coordinates(:,1)<500.1;
tmp2 = topology.coordinates(:,1)>499.9;
tmp3 = topology.coordinates(:,2)<500.1;
tmp4 = topology.coordinates(:,2)>499.9;
tmpNod = tmp1+tmp2+tmp3+tmp4;
vertNod = find(tmpNod == 4);
[vertNodZ,indNod] = sort(topology.coordinates(vertNod,3));


%find elemes in vertical symmetry axis
tmp1 = elems.cellCentroid(:,1)<450.1;
tmp2 = elems.cellCentroid(:,1)>449.9;
tmp3 = elems.cellCentroid(:,2)<550.1;
tmp4 = elems.cellCentroid(:,2)>449.9;
tmpEl = tmp1+tmp2+tmp3+tmp4;
vertEl = find(tmpEl == 4);
[vertElZ,indEl] = sort(elems.cellCentroid(vertEl,3));

timesInd = [2;3;4];
time_string = "Year  " + expTime(timesInd);
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.')

if isFEMBased(model,'Flow')
    pressPlot = expPress(vertNod(indNod),timesInd);
    figure(1)
    plot(pressPlot,vertNodZ)
    xlabel('Pressure [kPa]')
    ylabel('z (m)')
    legend(time_string)
elseif isFVTPFABased(model,'Flow')
    pressPlot = expPress(vertEl(indEl),timesInd);
    figure(1)
    plot(pressPlot,vertElZ)
    xlabel('Pressure [kPa]')
    ylabel('z (m)')
    legend(time_string)
end
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'DeepAcquifer_pressure', '.png');
exportgraphics(gcf,stmp,'Resolution',400)


dispPlot = expDispl(3*vertNod(indNod),timesInd);

figure(2)
plot(1000*dispPlot,vertNodZ);
xlabel('Vertical displacement (mm)')
ylabel('z (m)')
% xlim([0 50])
% ylim([-60 5])
legend(time_string)
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'DeepAcquifer_vertDisplacements', '.png');
exportgraphics(gcf,stmp,'Resolution',400)
