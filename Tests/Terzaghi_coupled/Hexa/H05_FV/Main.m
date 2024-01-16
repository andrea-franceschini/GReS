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
terzaghi_analytical;

%saving coordinates for later use
%save depht_H05.dat zvector  -ascii
%
% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
% Degree of freedom manager 
dofmanager = DoFManager(topology,model);

%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["dir_BC_flow_tetra.dat","dir_BCSurf_poro_tetra.dat","neuSurf_BC_poro_tetra.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid,dofmanager);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Some preprocessing stuff
%PreProc class has been removed in the last version
%indB is now defined in Elements Class
%getStiffMatrix is now a method of material SubClasses (Elastic, SSCM...)
%getDoFID is an external function in Discretizer repository
%pre = PreProc(grid,mat);
%
%reading vectors containing initial conditions for each node
% u0 = load('uIni05.dat');
% p0 = load('pIni05.dat');

% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
%file = 'initialconditions';
file = ["iniDisp.dat","iniPressure.dat"];
resState = State(model,grid,mat,file,GaussPts);

%manually assigning initial conditions before proper implementation
% resState.dispConv(3:3:end) = -u0fem'; %only DZ is fixed initially
% resState.dispCurr(3:3:end) = -u0fem'; %only DZ is fixed initially
% resState.pressure = p0fem;
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
plotObj1 = plot(pressplot,ptsY,'ko');
hold on
plotObj2 = plot(pfem(nodesP,:),ptsY,'k');
grid on
xlabel('Pressione (kPa)')
ylabel('z (m)')
legend([plotObj1(1),plotObj2(1)],{'Numerica','Analitica'});
title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')

figure(2)
plotObj1 = plot(-dispplot,topology.coordinates(nodesU,3),'ko');
hold on
plotObj2 = plot(ufem(nodesU,:),topology.coordinates(nodesU,3),'k');
grid on
xlabel('Spostamenti verticali (m)')
ylabel('z (m)')
title('h = 0.5 m \Delta t = 0.1 s \theta = 1.0')
legend([plotObj1(1),plotObj2(1)],{'Numerica','Analitica'});

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

%%












%
% % Compute the volume connected to each node
% volNod = zeros(topology.nNodes,1);
% if any(topology.cellVTKType == 12)
%   N1 = getBasisFinGPoints(elements.hexa);
% end
% for el=1:topology.nCells
%   top = topology.cells(el,1:topology.cellNumVerts(el));
%   if topology.cellVTKType(el) == 10 % Tetra
%     volNod(top) = volNod(top) + elems.vol(el)/topology.cellNumVerts(el);
%   elseif topology.cellVTKType(el) == 12 % Hexa
%     dJWeighed = getDerBasisFAndDet(elems.hexa,el,3);
%     volNod(top) = volNod(top)+ N1'*dJWeighed';
%   end
% end
% %
% % Edge length
% % ledge = zeros(topology.nCells,1);
% % for el = 1:topology.nCells
% %   comb = nchoosek(topology.cells(el,:),2);
% %   ledgeLoc = sqrt((topology.coordinates(comb(:,1),1)-topology.coordinates(comb(:,2),1)).^2 + ...
% %     (topology.coordinates(comb(:,1),2)-topology.coordinates(comb(:,2),2)).^2 + ...
% %     (topology.coordinates(comb(:,1),3)-topology.coordinates(comb(:,2),3)).^2);
% %   ledge(el) = max(ledgeLoc);
% % end
% %
% 
% % Analytical solution for flow problem
% %load('expData.mat');
% 
% qS = bound.getVals('neu_down', 1);
% qB = -qS(1);
% permMat = mat.getMaterial(2).getPermMatrix();
% kB = permMat(1,1);
% % fVec = bound.getVals('distrSource', 1);
% % fB = -fVec(1);
% fB = 0;
% pVec = bound.getVals('dir_top', 1);
% pB = pVec(1);
% len = max(topology.coordinates(:,3));
% pAnal = fB/(2*kB)*topology.coordinates(:,3).^2 + ...
%   qB/kB*topology.coordinates(:,3) + (pB-1/kB*((len^2)/2*fB+len*qB));
% errflow = (resState.pressure - pAnal).^2;
% errNormflow = sqrt(errflow'*volNod);
% 
% % Analytical solution_1D truss
% %load('expData.mat');
% fS = bound.getVals('neu_top', 1);
% fB = -fS(1);
% %permMat = mat.getMaterial(2).getPermMatrix();
% %kB = permMat(1,1);
% % fVec = bound.getVals('distrSource', 1);
% % fB = -fVec(1);
% E = mat.getMaterial(1).E;
% uVec = bound.getVals('dir_down', 1);
% uB = uVec(1);
% len = max(topology.coordinates(:,3));
% uAnal = fB/E*topology.coordinates(:,3);
% uz = resState.displ(3:3:end);
% errporo = (uz - uAnal).^2;
% errNormporo = sqrt(errporo'*volNod);
% 
% 
% delete(bound);