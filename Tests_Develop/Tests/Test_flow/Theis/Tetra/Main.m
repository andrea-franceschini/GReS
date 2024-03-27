close all;
clear;

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("SinglePhaseFlow_FEM");
%
% ----------------------- SIMULATION PARAMETERS ----------------------
fileName = "SimParam.dat";
simParam = SimulationParameters(model,fileName);
%
% ------------------------------  MESH -------------------------------
% Create the Mesh object
topology = Mesh();
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'Mesh_Theis_tetra.msh';
%
% Create an object of the Materials class and read the materials file
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
% Create an object of the "Elements" class and process the element properties
elems = Elements(topology);
%
% Wrap Mesh, Elements and Faces objects in a structure
grid = struct('topology',topology,'cells',elems,'faces',[]);
%
%------------------------ BOUNDARY CONDITIONS ------------------------
%
% Set the input file
% For transient simulations
fileName = ["dirNodeLatSurf_tetra.dat","neuNode.dat"];
% For steady-state simulations
% fileName = ["dirNodeLatSurf_tetra_SS.dat","neuNode.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);
%
%-------------------------- PREPROCESSING ----------------------------
%
% Set the "State" object. It contains all the vectors describing the state
% of the reservoir in terms of pressure, displacement, stress, ...
resState = State(model,grid,mat);
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
% Load the Theis model properties and numerical pressure solution
isSS = false;
if simParam.dtIni > simParam.tMax
  isSS = true;
end
load('../Analytical/TheisBenchSettings.mat');
load('expData.mat');
rDist = sqrt(topology.coordinates(:,1).^2 + topology.coordinates(:,2).^2);
if ~isSS
  [tVec,rVec] = meshgrid(expTime(2:end),rDist);
  u = rVec(:).^2*TModParam.S./(4*TModParam.T*tVec(:));
  W = expint(u);
  pAnal = TModParam.gamma*TModParam.q/(4*pi*TModParam.T)*W;
  pAnal = reshape(pAnal,length(rDist),length(expTime)-1);
  pAnal = TModParam.pIni-pAnal;
  err = (expPress(:,2:end) - pAnal).^2;
else
  expPress = expPress(:,1);
  pAnal = TModParam.pIni + (TModParam.q*TModParam.gamma)/ ...
    (2*pi*TModParam.T)*log(rDist/TModParam.R);
  err = (expPress - pAnal).^2;
end
%
err(bound.getDofs('nodeNeu'),:) = 0;
% remNod = (rDist <= 15);
% err(remNod,:) = 0;
errNorm = err'*volNod;
errNorm = sqrt(errNorm);
%
% Plots
nodPrint = [9464; 9394; 7216; 8520; 8202; 8196; 5834; 5836; 5906; 6034; ...
  6180; 8374; 8560; 7720; 8426; 5302; 5312; 5378; 5418; 5436; 5448; ...
  5470; 5602; 8664; 8662; 6264; 6430; 7008; 8766; 9238; 2588];
rNodPrint = sqrt(topology.coordinates(nodPrint,1).^2 + topology.coordinates(nodPrint,2).^2);
if isSS
  figure(1)
  fplot(@(r) TModParam.pIni + (TModParam.q*TModParam.gamma)/(2*pi*TModParam.T)* ...
    log(r/TModParam.R),[0.2,TModParam.R],'LineWidth',0.75,'DisplayName','Analytical Solution');
  hold on
  plot(rNodPrint,expPress(nodPrint),'LineStyle','none','Marker','square','DisplayName','Numerical Solution');
  title('Thiem well model');
  xlabel('r [m]');
  ylabel('p [kPa]');
  xlim([0 500]);
  ylim([300 400]);
  legend('Location','southeast');
%   set(legend,...
%     'Location','eastoutside');
%     'Position',[0.615595233803704 0.702936454484307 0.271071432862963 0.0830952398209343]);
% 'Location','eastoutside');
  hold off
  print -dpdf ValidateThiemSol.pdf;
else
  colorSeq = [     0    0.4470    0.7410;
                   0    0.4470    0.7410;
              0.8500    0.3250    0.0980;
              0.8500    0.3250    0.0980;
              0.9290    0.6940    0.1250;
              0.9290    0.6940    0.1250;
              0.4940    0.1840    0.5560;
              0.4940    0.1840    0.5560;
              0.4660    0.6740    0.1880;
              0.4660    0.6740    0.1880;
              0.3010    0.7450    0.9330;
              0.3010    0.7450    0.9330;
              0.6350    0.0780    0.1840;
              0.6350    0.0780    0.1840];
  colororder(colorSeq);
  figure(1)
  hold on
  for i=1:length(expTime)-1
    fplot(@(r) TModParam.pIni-TModParam.gamma*TModParam.q/ ...
      (4*pi*TModParam.T)*expint(r.^2*TModParam.S./(4*TModParam.T*expTime(i+1))), ...
      [0.2,TModParam.R],'LineWidth',0.75,'DisplayName',string(expTime(i+1)));
    plot(rNodPrint,expPress(nodPrint,i+1),'LineStyle','none','Marker','square', ...
      'DisplayName',string(expTime(i+1)));
  end
  title('Theis well model');
  xlabel('r [m]');
  ylabel('p [kPa]');
  xlim([0 500]);
  ylim([300 400]);
  legend
  legend('Location','southeast');
  print -dpdf ValidateTheisSol.pdf;
end
%
%
delete(bound);