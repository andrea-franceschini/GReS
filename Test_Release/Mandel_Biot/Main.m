close all;
clear;

%anal_path =  'C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Analytical_solution';

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
fileName = 'Mandel_H01.msh';
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
%
% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
dofmanager = DoFManager(topology,model);
%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["dirBotFixed.dat","dirLatFaceYPoro.dat","dirLatFaceXPoro.dat",...
    "topLoad.dat","dirFreeFaceFlow.dat"];
%
bound = Boundaries(fileName,model,grid,dofmanager);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
file = ["iniDisp.dat","iniPressure.dat"];
%
resState = State(model,grid,mat,file,GaussPts);
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
%Post processing using MAT-FILE
%list of nodes along vertical axis (with x,y=0)
tol = 0.001;
elemP1 = find(abs(elems.cellCentroid(:,3) - 0.05) < tol);
elemP2 = find(abs(elems.cellCentroid(:,2) - 0.025) < tol);
elemP = intersect(elemP1, elemP2);
nodesX1 = find(abs(topology.coordinates(:,2)-0.05)<tol) ;
nodesX2 = find(abs(topology.coordinates(:,3)-0.7)<tol);
nodesX = intersect(nodesX1,nodesX2);
nodesZ1 = find(abs(topology.coordinates(:,1)-0.6)<tol);
nodesZ2 = find(abs(topology.coordinates(:,2)-0.05)<tol);
nodesZ = intersect(nodesZ1,nodesZ2);
[coordsP,ind] = sort(elems.cellCentroid(elemP,1));
elemP = elemP(ind);
[coordsX,ind] = sort(topology.coordinates(nodesX,1));
nodesX = nodesX(ind);
[coordsZ,ind] = sort(topology.coordinates(nodesZ,3));
nodesZ = nodesZ(ind);

%Getting pressure and displacement solution for specified output times from MatFILE
press = printUtils.m.expPress;
disp = printUtils.m.expDispl;
pressNum = press(elemP,2:end);
dispXNum = disp(3*nodesX-2,2:end);
dispZNum = disp(3*nodesZ,2:end);



%Plotting solution
%Pressure
figure(1)
plotObj1 = plot(elems.cellCentroid(elemP,1),pressNum,'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(x,p,'k', 'LineWidth', 1);
xlabel('x (m)')
ylabel('Pressure (kPa)')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('C:\Users\Moretto\Documents\PHD\GReS\Reports\Presentation\Images\', 'Mandel_pressure', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

%Displacement DX
figure(2)
plotObj1 = plot(topology.coordinates(nodesX,1),dispXNum,'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(x,ux,'k', 'LineWidth', 1);
xlabel('x (m)')
ylabel('Displacement X (m)')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'southeast');
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('C:\Users\Moretto\Documents\PHD\GReS\Reports\Presentation\Images\', 'Mandel_UX', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

%Displacement DZ
figure(3)
plotObj1 = plot(dispZNum,topology.coordinates(nodesZ,3),'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(uz,z,'k', 'LineWidth', 1);
xlabel('Displacement Z (m)')
xlim([-10e-5 1e-5])
ylabel('z (m)')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});
grid on
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('C:\Users\Moretto\Documents\PHD\GReS\Reports\Presentation\Images\', 'Mandel_UZ', '.png');
exportgraphics(gcf,stmp,'Resolution',400)


%%
%Checking error norm 
% Compute the volume connected to each node
volNod = zeros(topology.nNodes,1);
if any(topology.cellVTKType == 12)
  N1 = getBasisFinGPoints(elems.hexa);
end
for el=1:topology.nCells
  top = topology.cells(el,1:topology.cellNumVerts(el));
  if topology.cellVTKType(el) == 10 % Tetra
    volNod(top) = volNod(top) + elems.vol(el)/topology.cellNumVerts(el);
  elseif topology.cellVTKType(el) == 12 % Hexa
    dJWeighed = getDerBasisFAndDet(elems.hexa,el,3);
    volNod(top) = volNod(top)+ N1'*dJWeighed';
  end
end


%errpress = sqrt(sum((analpress - press(:,2:end)).^2));
%normanal = sqrt(sum(analpress.^2));
%errRelpress = errpress./normanal;

%compute weighed error for the whole grid
errpress2 = (analpress - press(:,2:end)).^2;
errNormpress = sqrt(errpress2'*volNod);

