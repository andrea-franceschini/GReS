% close all;
% clear;

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
fileName = 'Test01_tetra.msh';
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

%xvector = topology.coordinates(:,1);
%zvector = topology.coordinates(:,3);

%mandel_analytical;
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
fileName = ["dir_Poro_fixed.dat","neuSurf_Poro_topLoad.dat","dir_Poro_Sym.dat","neuSurf_Flow_Impermeable.dat","dir_Flow_TopDrained.dat"];
%
% BClist = fileName;
% BCtable = {[1 2],5,4,3,[1,3,4],5};
% [entities_list,surf_list,loaded_surf] = writeBC(grid,BCtable,BClist);
%Create an object of the "Boundaries" class and read the boundary
%conditions
bound = Boundaries(fileName,model,grid);


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
NSolv = NonLinearSolver(model,simParam,grid,mat,bound,printUtils,resState);
%
% Solve the problem
[simState] = NSolv.NonLinearLoop();
%
% Finalize the print utility
printUtils.finalize()
%

%%
% -------------------------- POST PROCESSING ------------------------------

nodes = [143;215;287;359];
pressplot = expPress(nodes,:);
dispPlot = expDispl(nodes*3,:);
figure(1)
plot(expTime,pressplot,'k-o')
xlabel('Tempo (giorni)')
ylabel('Pressione (kPa)')
xlim([0 12])
ylim([-1 6])
grid on
figure(2)
plot(expTime,1000*dispPlot,'k-o')
xlabel('Tempo (giorni)')
ylabel('Spostamenti verticali (mm)')
xlim([0 12])
ylim([-40 5])
grid on

timesInd = [3 5 7 11 17];
nodes_subsidence = intersect(find(topology.coordinates(:,3) ==50),find(topology.coordinates(:,2) ==0));
x = sort(topology.coordinates(nodes_subsidence,1));
subs = sort(expDispl(3*nodes_subsidence,timesInd),'ascend');
set(0,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.')
figure(3)
plot(x,1000*subs);
xlabel('r (m)')
ylabel('Spostamento verticale (mm)')
xlim([0 50])
ylim([-60 5])
time_string = "Giorno  " + expTime(timesInd);
legend(time_string)

figure(4)
t = [0 1 6 7 100];
q = [0 -90000 -90000 0 0];
plot(t,q,'.-','LineWidth',1,'MarkerSize',15)
xlim([-1 10])
ylim([-100000 5000])
set(gca, 'YDir','reverse')
xlabel('Pressione (kPa)')
ylabel('Flusso prescritto (m^3/anno)')

figure(5)
t = [0 3 100];
q = [0 10 10];
plot(t,q,'.-','LineWidth',1,'MarkerSize',15)
xlim([-1 10])
ylim([-1 12])
xlabel('Tempo (giorni)')
ylabel('Carico agente (kPa)')

%grid on

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