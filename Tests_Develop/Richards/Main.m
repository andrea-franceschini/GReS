close all;
clear;

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("VariabSatFlow_FVTPFA");
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
fileName = 'Bench1D_hexa.msh';
%
% Import mesh data into the Mesh object
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
%
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology,GaussPts);
%
% Create an object of the "Faces" class and process the face properties
faces = Faces(model,topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
% fileName = ["neuSurfLeftFace_hexa2.dat","dirNodRightFace_hexa2.dat", ...
%   "volFBody_hexa2.dat"];
fileName = "dirBottom.dat";
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
dofmanager = DoFManager(topology,model);

bound = Boundaries(fileName,model,grid, dofmanager);
%linkBoundSurf2TPFAFace(model,bound,grid);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
fName = "iniPressure.dat";
resState = State(model,grid,mat,fName,GaussPts);
%
% Create and set the print utility
printUtils = OutState(model,mat,grid,'outTime.dat');
%
% Print the reservoir initial state
printUtils.printState(resState);
%
% ---------------------------- SOLUTION -------------------------------
%
% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam, dofmanager,grid,mat,bound, ...
  printUtils,resState,GaussPts);
%
% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%
% -------------------------- BENCHMARK ------------------------------
%
delete(bound);
%% -------------------------- BENCHMARK ------------------------------

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

tind = [9;21; 30; 37];
press = printUtils.m.expPress;
sw = printUtils.m.expSw;
t = printUtils.m.expTime;
t_max = t(end);
t = t(tind)/t_max;


tstr = strcat(num2str(t),' \cdotT');
%Getting pressure and saturation solution for specified time from MatFILE
pressplot = press(nodesP,tind);
swplot = sw(nodesP,tind);


%Plotting solution
if isFVTPFABased(model,'Flow')
    ptsY = elems.cellCentroid(nodesP,3);
else
    ptsY = topology.coordinates(nodesP,3);
end
figure(1)
plot(pressplot,ptsY,'.k-', 'LineWidth', 1, 'MarkerSize', 15);
hold on
xlabel('Pressure (kPa)')
ylabel('z(m)')
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 12)
% export figure with quality
stmp = strcat('C:\Users\Moretto\Documents\PHD\GReS\Reports\Presentation\Images\', 'Richards_pressure', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

figure(2)
plot(swplot,ptsY,'.k-', 'LineWidth', 1, 'MarkerSize', 15);
hold on
xlabel('Saturation')
ylabel('z(m)')
xlim([0 1])
str = strcat('t = ',tstr);
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 12)
% export figure with quality
stmp = strcat('C:\Users\Moretto\Documents\PHD\GReS\Reports\Presentation\Images\', 'Richards_staturation', '.png');
exportgraphics(gcf,stmp,'Resolution',400)
