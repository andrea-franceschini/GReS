close all;
clear;

%anal_path =  'C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Analytical_solution';

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType(["SinglePhaseFlow_FEM","Poromechanics_FEM"]);
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
fileName = 'Mandel_H01_tetra.msh';
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
%GaussPts = Gauss(12,2,3);
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology);

xvector = topology.coordinates(:,1);
zvector = topology.coordinates(:,3);

mandel_analytical;
%saving coordinates for later use
%save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Analytical_solution\xmesh.dat xvector  -ascii
%save C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Analytical_solution\zmesh.dat zvector  -ascii
%
% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["dirNodBotFacePoro.dat","dirNodLatFaceYPoro.dat","dirNodLatFaceXPoro.dat",...
    "neuSurfTopFacePoro.dat","neuNodTopBotLatFlow.dat","dirNodFreeFaceFlow.dat","neuSurfTopRigidYPoro.dat","neuSurfTopRigidXPoro.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);

BClist = fileName;
BCtable = {1;3;4;2;[2 3 4 1];5};
entities_list = writeBC(grid,BCtable,BClist);
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
%reading vectors containing initial conditions for each node, must change
%if the mesh change. The suffix represents the element dimension
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat);
%manually assigning initial conditions before proper implementation
resState.dispConv(3:3:end) = uz0fem'; 
resState.dispCurr(3:3:end) = uz0fem'; 
resState.dispConv(1:3:end) = ux0fem';
resState.dispConv(1:3:end) = ux0fem';
resState.pressure(1:end) = p0fem;
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
NSolv = NonLinearSolver(model,simParam,grid,mat,bound,printUtils,resState);
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
nodesP = find(topology.coordinates(:,2)+topology.coordinates(:,3)==0);
nodesX1 = find(abs(topology.coordinates(:,2)-0.5)<tol) ;
nodesX2 = find(abs(topology.coordinates(:,3)-0.5)<tol);
nodesX = intersect(nodesX1,nodesX2);
nodesZ1 = find(abs(topology.coordinates(:,1)-0.5)<tol);
nodesZ2 = find(abs(topology.coordinates(:,2)-0.5)<tol);
nodesZ = intersect(nodesZ1,nodesZ2);
[coordsP,ind] = sort(topology.coordinates(nodesP,1));
nodesP = nodesP(ind);
[coordsX,ind] = sort(topology.coordinates(nodesX,1));
nodesX = nodesX(ind);
[coordsZ,ind] = sort(topology.coordinates(nodesZ,3));
nodesZ = nodesZ(ind);

%Getting pressure and displacement solution for specified output times from MatFILE
press = printUtils.m.expPress;
disp = printUtils.m.expDispl;
pressplot = press(nodesP,2:end);
dispXplot = disp(3*nodesX-2,2:end);
dispZplot = disp(3*nodesZ,2:end);


%Getting analytical solution arrays for errors check
% analpress  = load('pAnal.dat');
% analDX = load('uxAnal.dat');
% analDZ = load('uzAnal.dat');

%getting position vectors for plots
% xAnal = load('xAnal.dat');
% zAnal = load('zAnal.dat');
%getting analytical solution arrays for plots
% analpressPlot  = analpress(nodesP,:);
% analDXPlot = analDX(nodesX,:);
% analDZPlot = analDZ(nodesZ,:);


%Plotting solution
%Pressure
figure(1)
plotObj1 = plot(topology.coordinates(nodesP,1),pressplot,'o');
hold on
plotObj2 = plot(topology.coordinates(nodesP,1),pfem(nodesP));
xlabel('x (m)')
ylabel('Pressure (kPa)')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});
title('h = 0.1 m \Delta t_{ini} = 0.01 s  \theta = 1.0')

%Displacement DX
figure(2)
plotObj1 = plot(topology.coordinates(nodesX,1),dispXplot,'o');
hold on
plotObj2 = plot(topology.coordinates(nodesX,1),uxfem(nodesX,1));
xlabel('X (m)')
ylabel('DX (m)')
title('h = 0.1 m \Delta t_{ini} = 0.01 s  \theta = 1.0')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});

%Displacement DZ
figure(3)
plotObj1 = plot(dispZplot,topology.coordinates(nodesZ,3),'o');
hold on
plotObj2 = plot(uzfem(nodesZ),topology.coordinates(nodesZ,3));
xlabel('Displacement Z (m)')
ylabel('Depht (m)')
title('h = 0.1 m \Delta t_{ini} = 0.01 s  \theta = 1.0')
legend([plotObj1(1),plotObj2(1)],{'Numerical','Analytical'});

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
errpress2 = (pfem - press(:,2:end)).^2;
errNormpress = sqrt(errpress2'*volNod);

errdispX2 = (uxfem - disp(1:3:end,2:end)).^2;
errNormDispX = sqrt(errdispX2'*volNod);

errdispZ2 = (uzfem - disp(3:3:end,2:end)).^2;
errNormDispZ = sqrt(errdispZ2'*volNod);





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