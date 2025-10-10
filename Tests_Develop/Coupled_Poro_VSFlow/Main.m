close all;
clear;

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType(["VariabSatFlow_FVTPFA","Poromechanics_FEM"]);
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
%terzaghi_analytical;

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
fileName = ["dir_BC_flow_tetra.dat","dir_BCSurf_poro_tetra.dat","neuSurf_BC_poro_tetra.dat"];
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


%Plotting solution
if isFVTPFABased(model,'Flow')
    ptsY = elems.cellCentroid(nodesP,3);
else
    ptsY = topology.coordinates(nodesP,3);
end
figure(1)
plotObj1 = plot(pressplot,ptsY,'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(pfem(nodesP,:),ptsY,'k', 'LineWidth', 1);
grid on
xlabel('Pressure (kPa)')
ylabel('z (m)')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'northeast');
%title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')
axis tight
xlim([0 10.2])

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 12)
% export figure with quality
stmp = strcat('C:\Users\Moretto\Documents\PHD\GReS\Reports\Presentation\Images', 'Terzaghi_pressure', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

figure(2)
plotObj1 = plot(-dispplot,topology.coordinates(nodesU,3),'k.', 'LineWidth', 1, 'MarkerSize', 15);
hold on
plotObj2 = plot(ufem(nodesU,:),topology.coordinates(nodesU,3),'k',  'LineWidth', 1);
grid on
xlabel('Vertical displacements (mm)')
ylabel('z (m)')
%title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'}, 'Location', 'southeast');

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif', 'FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif', 'FontSize', 12)
% export figure with quality
stmp = strcat('C:\Users\Moretto\Documents\PHD\GReS\Reports\Presentation\Images\', 'Terzaghi_disp', '.png');
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


% errpress = sqrt(sum((analpress - press(:,2:end)).^2));
% normanal = sqrt(sum(analpress.^2));
% errRelpress = errpress./normanal;

%pressure_error
errpress2 = (pfem - press(:,2:end)).^2;
errNormpressure = sqrt(errpress2'*elems.vol);

%displacement_error
errdisp2 = (ufem - disp(3:3:end,2:end)).^2;
errNormdisp = sqrt(errdisp2'*volNod);