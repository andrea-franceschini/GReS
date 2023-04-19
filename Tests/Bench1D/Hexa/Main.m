close all;
clear;

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("SinglePhaseFlow_FEM");
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "SimParam.dat";
simParam = SimulationParameters(fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
% Set the input file name
fileName = 'Bench1D_hexa1.msh';
%
% Import the mesh data into the Mesh object
topology.importGMSHmesh(fileName);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'materials.dat';
%
% Create an object of the Materials class and read the materials file
mat = Materials(fileName);
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
faces = Faces(topology,model);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',faces);
%
%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
fileName = ["dirElemRightFace_hexa1.dat","neuSurfLeftFace_hexa1.dat", ...
  "volFBody_hexa1.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName);
linkBoundSurf2TPFAFace(model,bound,grid);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Some preprocessing stuff
pre = PreProc(grid,mat,GaussPts);
%
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat,pre,GaussPts);
%
% Create and set the print utility
printUtils = OutState(model,grid,'outTime.dat');
%
% Print the reservoir initial state
printUtils.printState(resState);
%
% ---------------------------- SOLUTION -------------------------------
%
% Create the object handling the (nonlinear) solution of the problem
NSolv = NonLinearSolver(model,simParam,grid,mat,pre,bound, ...
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
load('expData.mat');
qS = bound.getVals('lFace', 1);
permMat = mat.getMaterial(2).getPermMatrix();
kB = permMat(1,1);
fVec = bound.getVals('distrSource', 1);
fB = -fVec(1);
pVec = bound.getVals('rFace', 1);
pB = pVec(1);
if model.isFEMBased('Flow')
  len = max(topology.coordinates(:,3));
  qB = -qS(1);
  pAnal = fB/(2*kB)*topology.coordinates(:,3).^2 + ...
    qB/kB*topology.coordinates(:,3) + (pB-1/kB*((len^2)/2*fB+len*qB));
  err = (expPress(:,1) - pAnal).^2;
  errNorm = sqrt(err'*volNod);
elseif model.isFVTPFABased('Flow')
  len = 10-0.125;
  qB = -sum(qS);
  pAnal = fB/(2*kB)*elems.cellCentroid(:,3).^2 + ...
    qB/kB*elems.cellCentroid(:,3) + (pB-1/kB*((len^2)/2*fB+len*qB));
  pAnal(end) = pB;
  err = (expPress(:,1) - pAnal).^2;
  errNorm = sqrt(err'*elems.vol);
end
%
delete(bound);