close all;
clear;

% -------------------------- SET THE PHYSICS -------------------------
model = ModelType("SinglePhaseFlow_FVTPFA");
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
fileName = 'Mesh_Theis_hexa_FVTPFA.msh';
%
% Import mesh data into the Mesh object
topology.importGMSHmesh(fileName);
% topology.cellTag = ones(length(topology.cellTag),1);
%
%----------------------------- MATERIALS -----------------------------
%
% Set the input file name
fileName = 'materialsList.dat';
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
%--------------------------- BOUNDARY CONDITIONS -------------------------
%
% fileName = ["dirNodeLatSurf_hexa_SS.dat","neuNode.dat"];
fileName = ["dirSurfLatSurf_hexa.dat","neuVol.dat"];
%
% Create an object of the "Boundaries" class and read the boundary
% conditions
bound = Boundaries(fileName,model,grid);
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
% Compute the volume connected to each node, if required
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
      volNod(top) = volNod(top) + N1'*dJWeighed';
    end
  end
  %
  rDist = sqrt(topology.coordinates(:,1).^2 + topology.coordinates(:,2).^2);
  itemPrint = [2650; 4968; 8914; 7752; 7310; 7312; 8686; 5828; 5822; 5820; ...
  9202; 7606; 7096; 6652; 6650; 8494; 8704; 8068; 5302; 5304; 6026; ...
  8446; 6924; 6804; 7978; 7944; 5644; 5646; 7976; 8598; 7404; 9360];
  rItemPrint = sqrt(topology.coordinates(itemPrint,1).^2 + ...
    topology.coordinates(itemPrint,2).^2);
elseif model.isFVTPFABased('Flow')
  rDist = sqrt(grid.cells.cellCentroid(:,1).^2 + grid.cells.cellCentroid(:,2).^2);
  itemPrint = [1592; 4157; 1172; 1685; 3719; 4436; 6305; 1235; 1220; ...
    3050; 3998; 1250; 767; 7466; 3572; 1862; 4100; 4229; 3812; 2387; ... %2387;
    3185; 4463; 3086; 3440; 3986; 3197; 4061; 89; 6680; 3227; 5069; ... 
    5318; 5273; 881; 512; 5354; 659; 8057];
  rItemPrint = sqrt(grid.cells.cellCentroid(itemPrint,1).^2 + ...
    grid.cells.cellCentroid(itemPrint,2).^2);
end
%
% Load the Theis model properties and numerical pressure solution
isSS = false;
if simParam.dtIni > simParam.tMax
  isSS = true;
end
load('../Analytical/TheisBenchSettings.mat');
load('expData.mat');
if ~isSS
  [tVec,rVec] = meshgrid(expTime(2:end),rDist);
  u = rVec(:).^2*TModParam.S./(4*TModParam.T*tVec(:));
  W = expint(u);
  pAnal = TModParam.gamma*TModParam.q/(4*pi*TModParam.T)*W;
  pAnal = reshape(pAnal,length(rDist),length(expTime) - 1);
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
if model.isFEMBased('Flow')
  errNorm = err'*volNod;
elseif model.isFVTPFABased('Flow')
  errNorm = err'*grid.cells.vol;
end
errNorm = sqrt(errNorm);
%
if isSS
  figure(1)
  fplot(@(r) TModParam.pIni + (TModParam.q*TModParam.gamma)/(2*pi*TModParam.T)* ...
    log(r/TModParam.R),[0.2,TModParam.R],'LineWidth',0.75,'DisplayName','Analytical Solution');
  hold on
  plot(rItemPrint,expPress(itemPrint),'LineStyle','none','Marker','square','DisplayName','Numerical Solution');
  if model.isFEMBased('Flow')
    title('Thiem well model - FEM');
  elseif model.isFVTPFABased('Flow')
    title('Thiem well model - FVTPFA');
  end
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
  if model.isFEMBased('Flow')
    print -dpdf ValidateThiemSol_FEM.pdf;
  elseif model.isFVTPFABased('Flow')
    print -dpdf ValidateThiemSol_FVTPFA.pdf;
  end
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
    plot(rItemPrint,expPress(itemPrint,i+1),'LineStyle','none','Marker','square', ...
      'DisplayName',string(expTime(i+1)));
  end
  if model.isFEMBased('Flow')
    title('Theis well model - FEM');
  elseif model.isFVTPFABased('Flow')
    title('Theis well model - FVTPFA');
  end
  xlabel('r [m]');
  ylabel('p [kPa]');
  xlim([0 500]);
  ylim([300 400]);
  legend
  legend('Location','southeast');
  if model.isFEMBased('Flow')
    print -dpdf ValidateTheisSol_FEM.pdf;
  elseif model.isFVTPFABased('Flow')
    print -dpdf ValidateTheisSol_FVTPFA.pdf;
  end
end
%
%
delete(bound);