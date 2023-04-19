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
fileName = 'Bench1D_tetra1.msh';
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
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology);
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
fileName = ["dirRightFace_tetra1.dat","neuSurfLeftFace_tetra1.dat", ...
  "volFBody_tetra1.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Some preprocessing stuff
pre = PreProc(grid,mat);
%
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat,pre);
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
NSolv = NonLinearSolver(model,simParam,grid,mat,pre,bound,printUtils,resState);
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
%
% Edge length
% ledge = zeros(topology.nCells,1);
% for el = 1:topology.nCells
%   comb = nchoosek(topology.cells(el,:),2);
%   ledgeLoc = sqrt((topology.coordinates(comb(:,1),1)-topology.coordinates(comb(:,2),1)).^2 + ...
%     (topology.coordinates(comb(:,1),2)-topology.coordinates(comb(:,2),2)).^2 + ...
%     (topology.coordinates(comb(:,1),3)-topology.coordinates(comb(:,2),3)).^2);
%   ledge(el) = max(ledgeLoc);
% end
%
% Analytical solution
load('expData.mat');
qS = bound.getVals('lFace', 1);
qB = -qS(1);
permMat = mat.getMaterial(2).getPermMatrix();
kB = permMat(1,1);
fVec = bound.getVals('distrSource', 1);
fB = -fVec(1);
% fB = 0;
pVec = bound.getVals('rFace', 1);
pB = pVec(1);
len = max(topology.coordinates(:,3));
pAnal = fB/(2*kB)*topology.coordinates(:,3).^2 + ...
  qB/kB*topology.coordinates(:,3) + (pB-1/kB*((len^2)/2*fB+len*qB));
err = (expPress(:,1) - pAnal).^2;
errNorm = sqrt(err'*volNod);
%
delete(bound);