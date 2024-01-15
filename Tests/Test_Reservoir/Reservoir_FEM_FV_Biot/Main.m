close all;
clear;

%adding path
rmpath(genpath('C:\Users\Moretto\Documents\PHD\GReS\GReS'))
addpath(genpath('C:\Users\Moretto\Documents\PHD\GReS\GReS\Code'));
addpath(genpath(pwd));
%anal_path =  'C:\Users\Moretto\Documents\UNIPD\Tesi_magistrale\Code_18_07\GReS\Tests\Mandel_coupled\Analytical_solution';
tic
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
fileName = 'ReservoirTest_Hexa.msh';
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'materialsListElastic.dat';
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
% Create an object of the "Faces" class and process the face properties
faces = Faces(model, topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
%----------------------------- DOF MANAGER -----------------------------
fileName = 'dof.dat';
dofmanager = DoFManager(topology, model, fileName);

%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["bottom_fixed.dat","flux.dat","lateral_fix.dat","Impermeable.dat"];
%
% BClist = fileName;
% BCtable = {1,2,[1 3 4],3,4};
% [entities_list,surf_list] = writeBC(grid,BCtable,BClist);
%Create an object of the "Boundaries" class and read the boundary
%conditions
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
%reading vectors containing initial conditions for each node, must change
%if the mesh change. The suffix represents the element dimension
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat,GaussPts);
%manually assigning initial conditions before proper implementation
% resState.dispConv(3:3:end) = uz0fem'; 
% resState.dispCurr(3:3:end) = uz0fem'; 
% resState.dispConv(1:3:end) = ux0fem';
% resState.dispConv(1:3:end) = ux0fem';
% resState.pressure(1:end) = p0fem;
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
time = toc;

%% POST PROCESSING
%loading output file (inside mandel directory for some reason)
load('expData.mat')

%find nodes in vertical symmetry axis
tmp1=topology.coordinates(:,1)<500.1;
tmp2 = topology.coordinates(:,1)>499.9;
tmp3 = topology.coordinates(:,2)<500.1;
tmp4 = topology.coordinates(:,2)>499.9;
tmpNod = tmp1+tmp2+tmp3+tmp4;
vertNod = find(tmpNod == 4);
[vertNodZ,indNod] = sort(topology.coordinates(vertNod,3));


%find elemes in vertical symmetry axis
tmp1 = elems.cellCentroid(:,1)<450.1;
tmp2 = elems.cellCentroid(:,1)>449.9;
tmp3 = elems.cellCentroid(:,2)<550.1;
tmp4 = elems.cellCentroid(:,2)>449.9;
tmpEl = tmp1+tmp2+tmp3+tmp4;
vertEl = find(tmpEl == 4);
[vertElZ,indEl] = sort(elems.cellCentroid(vertEl,3));

timesInd = [2;5;8];
time_string = "Year  " + expTime(timesInd);
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.')

if isFEMBased(model,'Flow')
    pressPlot = expPress(vertNod(indNod),timesInd);
    figure(1)
    plot(pressPlot,vertNodZ)
    xlabel('Pressione [kPa]')
    ylabel('z (m)')
    legend(time_string)
elseif isFVTPFABased(model,'Flow')
    pressPlot = expPress(vertEl(indEl),timesInd);
    figure(1)
    plot(pressPlot,vertElZ)
    xlabel('Pressione [kPa]')
    ylabel('z (m)')
    legend(time_string)
end
dispPlot = expDispl(3*vertNod(indNod),timesInd);

figure(2)
plot(1000*dispPlot,vertNodZ);
xlabel('Spostamento verticale (mm)')
ylabel('z (m)')
% xlim([0 50])
% ylim([-60 5])
legend(time_string)



%%
%Checking error norm 
% Compute the volume connected to each node
% volNod = zeros(topology.nNodes,1);
% if any(topology.cellVTKType == 12)
%   N1 = getBasisFinGPoints(elems.hexa);
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


%errpress = sqrt(sum((analpress - press(:,2:end)).^2));
%normanal = sqrt(sum(analpress.^2));
%errRelpress = errpress./normanal;

%compute weighed error for the whole grid
% errpress2 = (pfem - press(:,2:end)).^2;
% errNormpress = sqrt(errpress2'*volNod);
% 
% errdispX2 = (uxfem - disp(1:3:end,2:end)).^2;
% errNormDispX = sqrt(errdispX2'*volNod);
% 
% errdispZ2 = (uzfem - disp(3:3:end,2:end)).^2;
% errNormDispZ = sqrt(errdispZ2'*volNod);

%diagram to find stationary flow
% t=[10 50 100 500 1000];
% p=[678 2580 4239 11852 17741.3];
% plot(t,p)



%%
%
% % Compute the volume connected to each node
volNod = zeros(topology.nNodes,1);
if any(topology.cellVTKType == 12)
  N1 = getBasisFinGPoints(elements.hexa);
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