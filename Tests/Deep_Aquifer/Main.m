close all;
% clear;

profile on

scriptFullPath = mfilename('fullpath');
scriptDir = fileparts(scriptFullPath);
cd(scriptDir);
warning('off','MATLAB:nearlySingularMatrix');

% <<<<<<< HEAD
% shortcut to define a model using a unique xml file
% useful when dealing with many domains
domain = buildModel('domain.xml');

% perform a fully coupled simulation
solver = FCSolver(domain);
[simState] = solver.NonLinearLoop();

% Finalize the print utility
domain.outstate.finalize()

% =======
% % List the physical models activated in the simulation and their
% % discretization scheme
% model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);
% % model = ModelType(["SinglePhaseFlow_FEM","Poromechanics_FEM"]);
% 
% % Create object containing simulation parameters
% fileName = "simParam.dat";
% simParam = SimulationParameters(fileName,model);
% 
% % Create the Mesh object
% topology = Mesh();
% 
% % Set the input file name
% fileName = 'ReservoirTest_Hexa.msh';
% 
% % Import the mesh data into the Mesh object
% topology.importGMSHmesh(fileName);
% 
% % Set the material input file name
% fileName = 'materialsListElastic.dat';
% %
% % Create an object of the Materials class and read the materials file
% mat = Materials(model,fileName);
% 
% % Define Gauss points
% GaussPts = Gauss(12,2,3);
% 
% % Create an object of the "Elements" class and process the element properties
% elems = Elements(topology,GaussPts);
% 
% % Create an object of the "Faces" class and process the face properties
% faces = Faces(model, topology);
% %
% % Wrap Mesh, Elements and Faces objects in a structure
% grid = struct('topology',topology,'cells',elems,'faces',faces);
% %
% % create a DoFManager object
% DoFfileName = 'dof.dat';
% dofmanager = DoFManager(topology, model, DoFfileName);
% 
% % Create object handling construction of Jacobian and rhs of the model
% linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);
% 
% % Build a structure storing variable fields at each time step
% state = linSyst.setState();
% 
% % Create and set the print utility
% printUtils = OutState(model,topology,'outTime.dat','folderName','Output_DeepAquifer');
% 
% % Write BC files programmatically with function utility 
% %fileName = setDeepAquiferBC('BCs',time,flux,topology,DoFfileName);
% if ~isempty(DoFfileName)
%     fileName = ["BCs/bottom_fixed.dat","BCs/flux.dat","BCs/lateral_fix.dat","BCs/Impermeable.dat"];
% else
%     fileName = ["BCs/bottom_fixed.dat","BCs/flux.dat","BCs/lateral_fix.dat","BCs/Impermeable.dat", "BCs/top_drained.dat"];
% end
% 
% bound = Boundaries(fileName,model,grid);
% 
% % perform a fully coupled simulation
% solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);
% [simState] = solver.NonLinearLoop();
% 
% % Finalize the print utility
% printUtils.finalize()
% >>>>>>> 1dfffa00097f21a2e1d34699913ab58ea5431391

%% POST PROCESSING

image_dir = strcat(pwd,'/Images');
if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir Images
else
    mkdir Images
end

% <<<<<<< HEAD
expPress = domain.outstate.results.expPress;
expDispl = domain.outstate.results.expDispl;
expTime = domain.outstate.results.expTime;

topol = domain.grid.topology;

%find nodes in vertical symmetry axis
tmp1=topol.coordinates(:,1)<500.1;
tmp2 = topol.coordinates(:,1)>499.9;
tmp3 = topol.coordinates(:,2)<500.1;
tmp4 = topol.coordinates(:,2)>499.9;
tmpNod = tmp1+tmp2+tmp3+tmp4;
vertNod = find(tmpNod == 4);
[vertNodZ,indNod] = sort(topol.coordinates(vertNod,3));


%find elemes in vertical symmetry axis
tmp1 = topol.cellCentroid(:,1)<450.1;
tmp2 = topol.cellCentroid(:,1)>449.9;
tmp3 = topol.cellCentroid(:,2)<550.1;
tmp4 = topol.cellCentroid(:,2)>449.9;
tmpEl = tmp1+tmp2+tmp3+tmp4;
vertEl = find(tmpEl == 4);
[vertElZ,indEl] = sort(topol.cellCentroid(vertEl,3));
% =======
% % expPress = printUtils.results.expPress;
% % expDispl = printUtils.results.expDispl;
% % expTime = printUtils.results.expTime;
% 
% % Small modification - for the growning grid
% nrep = length(printUtils.results(:,:));
% nvars = length(printUtils.results(2,:).expPress);
% expPress = zeros(nvars,nrep);
% nvars = length(printUtils.results(2,:).expDispl);
% expDispl = zeros(nvars,nrep);
% nvars = length(printUtils.results(2,:).expTime);
% expTime = zeros(nvars,nrep);
% for i=2:nrep
%    expPress(:,i) = printUtils.results(i,:).expPress;
%    expDispl(:,i) = printUtils.results(i,:).expDispl;
%    expTime(:,i) = printUtils.results(i,:).expTime;
% end
% 
% 
% 
% %find nodes in vertical symmetry axis
% tmp1=topology.coordinates(:,1)<500.1;
% tmp2 = topology.coordinates(:,1)>499.9;
% tmp3 = topology.coordinates(:,2)<500.1;
% tmp4 = topology.coordinates(:,2)>499.9;
% tmpNod = tmp1+tmp2+tmp3+tmp4;
% vertNod = find(tmpNod == 4);
% [vertNodZ,indNod] = sort(topology.coordinates(vertNod,3));
% 
% 
% %find elemes in vertical symmetry axis
% tmp1 = elems.cellCentroid(:,1)<450.1;
% tmp2 = elems.cellCentroid(:,1)>449.9;
% tmp3 = elems.cellCentroid(:,2)<550.1;
% tmp4 = elems.cellCentroid(:,2)>449.9;
% tmpEl = tmp1+tmp2+tmp3+tmp4;
% vertEl = find(tmpEl == 4);
% [vertElZ,indEl] = sort(elems.cellCentroid(vertEl,3));
% >>>>>>> 1dfffa00097f21a2e1d34699913ab58ea5431391

timesInd = [2;3;4];
time_string = "Year  " + expTime(timesInd);
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.')

% <<<<<<< HEAD
if isFEMBased(domain.model,'Flow')
% =======
% if isFEMBased(model,'Flow')
% >>>>>>> 1dfffa00097f21a2e1d34699913ab58ea5431391
    pressPlot = expPress(vertNod(indNod),timesInd);
    figure(1)
    plot(pressPlot,vertNodZ)
    xlabel('Pressure [kPa]')
    ylabel('z (m)')
    legend(time_string)
% <<<<<<< HEAD
elseif isFVTPFABased(domain.model,'Flow')
% =======
% elseif isFVTPFABased(model,'Flow')
% >>>>>>> 1dfffa00097f21a2e1d34699913ab58ea5431391
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
