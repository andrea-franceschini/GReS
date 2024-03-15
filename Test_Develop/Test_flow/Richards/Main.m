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
fileName = 'Bench1D_hexa3.msh';
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
fileName = "dirBottom3.dat";
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
resState = State(model,grid,mat,GaussPts);
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

tind = [3;5;10;20;30;37];
press = printUtils.m.expPress;
sw = printUtils.m.expSw;
t = printUtils.m.expTime;
t = t(tind);
tstr = num2str(t);
%Getting pressure and displacement solution for specified time from MatFILE
pressplot = press(nodesP,tind);
swplot = sw(nodesP,tind);


%Plotting solution
if isFVTPFABased(model,'Flow')
    ptsY = elems.cellCentroid(nodesP,3);
else
    ptsY = topology.coordinates(nodesP,3);
end
figure(1)
plot(pressplot,ptsY,'-o');
hold on
xlabel('Pressione (kPa)')
ylabel('z (m)')
legend(tstr)


figure(2)
plot(-swplot,ptsY,'-o');
hold on
xlabel('Saturazione')
ylabel('z (m)')
str = strcat('t = ',tstr);
legend(str)
