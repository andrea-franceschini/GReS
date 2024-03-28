close all;
clear;

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
fileName = 'TerzaghiH05_hexa.msh';
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
Terzaghi_analytical(topology, mat, 10)
%------------------------------ ELEMENTS -----------------------------
%
GaussPts = Gauss(12,2,3);
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);

% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
% Degree of freedom manager 
fname = 'dof.dat';
dofmanager = DoFManager(topology,model,fname);

%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["dir_BC_flow.dat","dir_BCSurf_poro.dat","neuSurf_BC_poro.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid,dofmanager);

%file = 'initialconditions';
if isFEMBased(model,'Flow')
    file = ["iniDisp.dat","iniPressureFEM.dat"];
else
    file = ["iniDisp.dat","iniPressure.dat"];
end

resState = State(model,grid,mat,file,GaussPts);

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
press = printUtils.m.expPress;
disp = printUtils.m.expDispl;
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
