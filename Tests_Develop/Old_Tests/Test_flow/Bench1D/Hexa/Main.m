close all;
clear;

%% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("SinglePhaseFlow_FVTPFA");
%
%% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "SimParam.dat";
simParam = SimulationParameters(model,fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'Bench1D_hexa4.msh';
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
%
% For FEM
% fileName = ["neuSurfLeftFace_hexa1.dat","dirNodRightFace_hexa1.dat", ...
%   "volFBody_hexa1.dat"];
%
% For FVTPFA
fileName = ["neuSurfLeftFace_hexa4.dat","dirSurfRightFace_hexa4.dat", ...
  "volFBody_hexa4.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
<<<<<<< HEAD:Tests/Test_flow/Bench1D/Hexa/Main.m
bound = Boundaries(fileName,model,grid);
linkBoundSurf2TPFAFace(model,bound,grid);
=======
bound = Boundaries(fileName);
%linkBoundSurf2TPFAFace(model,bound,grid);
>>>>>>> 2fd0824256e1ec6d23d32b48da985c601d732b96:Tests/Bench1D/Hexa/Main.m
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
NSolv = NonLinearSolver(model,simParam,grid,mat,bound, ...
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
% Compute the volume connected to each node
if model.isFEMBased('Flow')
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
end
%
% Benchmark
% if 5<1
load('expData.mat');
qS = bound.getVals('lFace', 1);
permMat = mat.getMaterial(1).PorousRock.getPermMatrix();
kB = permMat(1,1);
fVec = bound.getVals('distrSource', 1);
fB = -fVec(1);
pVec = bound.getVals('rFace', 1);
pB = pVec(1);
len = max(topology.coordinates(:,3));
qB = -qS(1);
pAnal = @(x) fB/(2*kB)*(x.^2) + qB/kB*x + (pB-1/kB*((len^2)/2*fB+len*qB));
% intPAnal = @(x) fB/(6*kB)*(x.^3) + qB/(2*kB).*(x^2) + (pB-1/kB*((len^2)/2*fB+len*qB)).*x;
if model.isFEMBased('Flow')
  pAnalNod = pAnal(topology.coordinates(:,3));
  err = (expPress(:,1) - pAnalNod).^2;
  errNorm = sqrt(err'*volNod);
elseif model.isFVTPFABased('Flow')
%   l = 0.125;
%   len = 10-l;
%   qB = -qS'*elems.vol(1:length(qS));
  %%%%
  %%%
%   qB = qS(1);
%   fB = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   errNorm = 0;
%   pAnalElem = zeros(grid.topology.nCells,1);
%   for el=1:grid.topology.nCells
%     top = grid.topology.cells(el,1:grid.topology.cellNumVerts(el));
%     zMax = max(grid.topology.coordinates(top,3));
%     zMin = min(grid.topology.coordinates(top,3));
%     pAnalElem(el) = (intPAnal(zMax)-intPAnal(zMin))/(zMax-zMin);
%     diffInt = pAnalElem(el) - expPress(el,1);
%     errNorm = errNorm + diffInt^2*elems.vol(el);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pAnalElem = pAnal(elems.cellCentroid(:,3));
    err = (expPress(:,1) - pAnalElem).^2;
    errNorm = sqrt(err'*elems.vol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   end
% errNorm = sqrt(errNorm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% end
%
delete(bound);