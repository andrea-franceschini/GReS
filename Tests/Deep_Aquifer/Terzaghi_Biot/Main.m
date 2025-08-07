close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

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
fileName = 'Mesh/Column.msh';
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
GaussPts = Gauss(12,2,3);
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

%saving z_vector coordinates for analytical solution calculation
nodez = topology.coordinates(:,3);
if isFVTPFABased(model,'Flow')
    cellz = elems.cellCentroid(:,3);
else
    cellz = nodez;
end

%calling analytical solution script
Terzaghi_analytical(topology, mat, 10)

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
% Degree of freedom manager 
%fname = 'dof.dat';
dofmanager = DoFManager(topology,model);

% Create object handling construction of Jacobian and rhs of the model
linSyst = Discretizer(model,simParam,dofmanager,grid,mat,GaussPts);

% Build a structure storing variable fields at each time step
state = linSyst.setState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','Output_Terzaghi');


%------------------------ BOUNDARY CONDITIONS ------------------------
% Write BC files programmatically with function utility 
F = -10; % vertical force
setTerzaghiBC('BCs',F,topology);

% Collect BC input file in a list
fileName = ["BCs/dirFlowTop.dat","BCs/newPorotop.dat",...
   "BCs/dirPoroLatY.dat","BCs/dirPoroLatX.dat","BCs/dirPoroBottom.dat"];
%
% Create an object of the "Boundaries" class 
bound = Boundaries(fileName,model,grid);

% In this version of the code, the user can assign initial conditions only
% manually, by directly modifying the entries of the state structure. 
% In this example, we use a user defined function to apply Terzaghi initial
% conditions to the state structure
state = applyTerzaghiIC(state,mat,topology,F);

% Print model initial state
printUtils.printState_new(linSyst,state);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a function
Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,state,linSyst,GaussPts);
%
% Solve the problem
[simState] = Solver.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%
%%
% -------------------------- BENCHMARK ------------------------------
image_dir = strcat(pwd,'/Images');
if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir Images
else
    mkdir Images
end

load("Terzaghi_Analytical.mat");
%Post processing using MAT-FILE 
%nodes vector contain list of nodes along vertical axis (with x,y=0)
nodesU = find(topology.coordinates(:,1)+topology.coordinates(:,2)==0);
[~,ind] = sort(topology.coordinates(nodesU,3));
nodesU = nodesU(ind);


% elem vector containing elements centroid along vertical axis
if isFEMBased(model,'Flow')
    nodesP = nodesU;
else
    nodesP = find(elems.cellCentroid(:,1) + elems.cellCentroid(:,2) < 0.51);
    [~,ind] = sort(elems.cellCentroid(nodesP,3));
    nodesP = nodesP(ind);
end

%Getting pressure and displacement solution for specified time from MatFILE
press = printUtils.results.expPress;
disp = printUtils.results.expDispl;
pressplot = press(nodesP,2:end);
dispplot = disp(3*nodesU,2:end);

H = max(topology.coordinates(:,3));
p0 = max(press(:,1));


%Plotting solution
if isFVTPFABased(model,'Flow')
    ptsY = elems.cellCentroid(nodesP,3);
else
    ptsY = topology.coordinates(nodesP,3);
end
figure(1)
plotObj1 = plot(pressplot/p0,ptsY/H,'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(p/p0,z/H,'k-', 'LineWidth', 1);
grid on
xlabel('p/p_0')
ylabel('z/H')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'northeast');
%title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')
axis tight
xlim([0 1.02])
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'Terzaghi_pressure', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

figure(2)
plotObj1 = plot(-dispplot/H,topology.coordinates(nodesU,3)/H,'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(u/H,z/H,'k-',  'LineWidth', 1);
grid on
xlabel('u_z/H')
ylabel('z/H')
%title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'southeast');

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'Terzaghi_disp', '.png');
exportgraphics(gcf,stmp,'Resolution',400)