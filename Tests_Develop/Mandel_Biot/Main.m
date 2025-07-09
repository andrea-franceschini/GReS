close all;
clear;

% Get the full path of the currently executing file
scriptFullPath = mfilename('fullpath');

% Extract the directory containing the script
scriptDir = fileparts(scriptFullPath);

% Change the current directory to the script's directory
cd(scriptDir);

% Set physical models 
model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);


% Set parameters of the simulation
fileName = "simParam.dat";
simParam = SimulationParameters(fileName,model);

% Create the Mesh object
topology = Mesh();

% Set the mesh input file name
fileName = 'Mandel_Mesh.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);

% Create an object of the Materials class and read the materials file
fileName = 'materialsList.dat';
mat = Materials(model,fileName);

% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,2);

%calling analytical solution script
Mandel_Analytical(topology, mat, 10)

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
linSyst = Discretizer(model,simParam,dofmanager,grid,mat);

% Build a structure storing variable fields at each time step
linSyst.initState();

% Create and set the print utility
printUtils = OutState(model,topology,'outTime.dat','folderName','Output_Mandel','flagMatFile',true);


% Write BC files programmatically with function utility 
F = -10; % vertical force
bcList = setMandelBC('BCs',F,topology);

% Create an object of the "Boundaries" class 
bound = Boundaries(bcList,model,grid);

% In this version of the code, the user can assign initial conditions only
% manually, by directly modifying the entries of the state structure. 
% In this example, we use a user defined function to apply Terzaghi initial
% conditions to the state structure
state = applyMandelIC(linSyst.state,mat,topology,F);

% Print model initial state
printState(printUtils,linSyst);

% The modular structure of the discretizer class allow the user to easily
% customize the solution scheme. 
% Here, a built-in fully implict solution scheme is adopted with class
% FCSolver. This could be simply be replaced by a user defined function
Solver = FCSolver(model,simParam,dofmanager,grid,mat,bound,printUtils,linSyst);
%
% Solve the problem
[simState] = Solver.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%
%% POST PROCESSING

image_dir = strcat(pwd,'/Images');
if isfolder(image_dir)
    rmdir(image_dir,"s")
    mkdir Images
else
    mkdir Images
end

%Post processing using MAT-FILE
%list of nodes along vertical axis (with x,y=0)
tol = 0.001;
elemP1 = find(abs(topology.cellCentroid(:,3) - 0.05) < tol);
elemP2 = find(abs(topology.cellCentroid(:,2) - 0.025) < tol);
elemP = intersect(elemP1, elemP2);
nodesX1 = find(abs(topology.coordinates(:,2)-0.05)<tol) ;
nodesX2 = find(abs(topology.coordinates(:,3)-0.7)<tol);
nodesX = intersect(nodesX1,nodesX2);
nodesZ1 = find(abs(topology.coordinates(:,1)-0.6)<tol);
nodesZ2 = find(abs(topology.coordinates(:,2)-0.05)<tol);
nodesZ = intersect(nodesZ1,nodesZ2);
[coordsP,ind] = sort(topology.cellCentroid(elemP,1));
elemP = elemP(ind);
[coordsX,ind] = sort(topology.coordinates(nodesX,1));
nodesX = nodesX(ind);
[coordsZ,ind] = sort(topology.coordinates(nodesZ,3));
nodesZ = nodesZ(ind);

%Getting pressure and displacement solution for specified output times from MatFILE
press = printUtils.results.expPress;
disp = printUtils.results.expDispl;
pressNum = press(elemP,2:end);
dispXNum = disp(3*nodesX-2,2:end);
dispZNum = disp(3*nodesZ,2:end);

% Loading analytical solution MAT-FILE into workspace
load("Mandel_Analytical.mat");

%Plotting solution
%Pressure
figure(1)
plotObj1 = plot(topology.cellCentroid(elemP,1),pressNum,'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(x,p,'k-', 'LineWidth', 1);
xlabel('x (m)')
ylabel('Pressure (kPa)')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'Mandel_pressure', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

%Displacement DX
figure(2)
plotObj1 = plot(topology.coordinates(nodesX,1),dispXNum,'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(x,ux,'k-', 'LineWidth', 1);
xlabel('x (m)')
ylabel('Displacement u_x (m)')
ylim([0 7.e-5])
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'southeast');
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'Mandel_UX', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

%Displacement DZ
figure(3)
plotObj1 = plot(dispZNum,topology.coordinates(nodesZ,3),'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(uz,z,'k-', 'LineWidth', 1);
xlabel('Displacement u_z (m)')
xlim([-10e-5 1e-5])
ylabel('z (m)')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'Mandel_UZ', '.png');
exportgraphics(gcf,stmp,'Resolution',400)
