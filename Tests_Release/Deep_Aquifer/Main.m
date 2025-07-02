close all;
clear;

profile on
scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);
warning('off','MATLAB:nearlySingularMatrix');

% List the physical models activated in the simulation and their
% discretization scheme
model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);

% Create object containing simulation parameters
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);

% Create the Mesh object
topology = Mesh();

% Set the input file name
fileName = 'ReservoirTest_Hexa.msh';

% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

% Set the material input file name
fileName = 'materialsListElastic.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(model,fileName);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
% create a DoFManager object
DoFfileName = 'dof.dat';
dofmanager = DoFManager(topology, model, DoFfileName);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat);

% Build a structure storing variable fields at each time step
linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','Output_DeepAquifer','flagMatFile',true);

% Write BC files programmatically with function utility 
%fileName = setDeepAquiferBC('BCs',time,flux,topology,DoFfileName);
if ~isempty(DoFfileName)
    fileName = ["BCs/bottom_fixed.dat","BCs/flux.dat","BCs/lateral_fix.dat","BCs/Impermeable.dat"];
else
    fileName = ["BCs/bottom_fixed.dat","BCs/flux.dat","BCs/lateral_fix.dat","BCs/Impermeable.dat", "BCs/top_drained.dat"];
end

bound = Boundaries(fileName,model,grid);

% perform a fully coupled simulation
solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,linSyst);
[simState] = solver.NonLinearLoop();

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

expPress = printUtils.results.expPress;
expDispl = printUtils.results.expDispl;
expTime = printUtils.results.expTime;

%find nodes in vertical symmetry axis
tmp1=topology.coordinates(:,1)<500.1;
tmp2 = topology.coordinates(:,1)>499.9;
tmp3 = topology.coordinates(:,2)<500.1;
tmp4 = topology.coordinates(:,2)>499.9;
tmpNod = tmp1+tmp2+tmp3+tmp4;
vertNod = find(tmpNod == 4);
[vertNodZ,indNod] = sort(topology.coordinates(vertNod,3));


%find elemes in vertical symmetry axis
tmp1 = topology.cellCentroid(:,1)<450.1;
tmp2 = topology.cellCentroid(:,1)>449.9;
tmp3 = topology.cellCentroid(:,2)<550.1;
tmp4 = topology.cellCentroid(:,2)>449.9;
tmpEl = tmp1+tmp2+tmp3+tmp4;
vertEl = find(tmpEl == 4);
[vertElZ,indEl] = sort(topology.cellCentroid(vertEl,3));

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
