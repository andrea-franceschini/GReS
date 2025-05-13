%clear; close all; 

%% MacroElement approach to spot violation of inf-sup condition

% Input: specify number of elements for top (master) and bottom (slave)
% side of the mesh.

% minimum number of elements is 2 for both sides

% all boundary nodes have fixed displacements.
% each element has fixed size

nG = 6;
E = 1; 
nu = 0;
integration = "SB"; % RBF,SB,P0
mult_type = 'dual';
bound = true;

Dmat = zeros(3);
Dmat([1 5]) = 1-nu;
Dmat([2 4]) = nu;
Dmat(9) = 0.5*(1-2*nu);
Dmat = (E/((1+nu)*(1-2*nu)))*Dmat;

% Set the input file name
gaussQuad = Gauss(12,2,2);

NM0X = 2;
NM0Y = 2;
NS0X = 4;
NS0Y = 4;
nR = 4;
h = zeros(nR,1);
infSup = zeros(nR,1);

for iref = 1:nR
% read mesh
NMX = NM0X*2^(iref-1);
NSX = NS0X*2^(iref-1);
NMY = NM0Y*2^(iref-1);
NSY = NS0Y*2^(iref-1);
h(iref) = 1/NSX;
getPatchMesh('Mesh/master.geo','master',NMX,NMY);
getPatchMesh('Mesh/slave.geo','slave',NSX,NSY);
[slaveMesh,masterMesh] = deal(Mesh(),Mesh()); 
slaveMesh.importGMSHmesh('Mesh/slave.msh')
masterMesh.importGMSHmesh('Mesh/master.msh');
yMax = max(slaveMesh.coordinates(:,2));
masterMesh.coordinates(:,2) = masterMesh.coordinates(:,2) + yMax;


Em = E*ones(masterMesh.nSurfaces,1); 
Es = E*ones(slaveMesh.nSurfaces,1); 
% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh,gaussQuad);
elemsSlave = Elements(slaveMesh,gaussQuad);


KMaster = stiff(masterMesh, elemsMaster, Dmat, gaussQuad);
KSlave = stiff(slaveMesh, elemsSlave, Dmat, gaussQuad);


mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
dof = DofMap(masterMesh,mortar.nodesMaster,slaveMesh,mortar.nodesSlave);


% compute mortar operator1
[Drbf, Mrbf] = mortar.computeMortarRBF(nG,4,'gauss','standard');
[Dsb, Msb] = mortar.computeMortarSegmentBased(nG,mult_type);
[Dconst,Mconst,Dconsistent] = mortar.computeMortarConstant(nG,4);
%
dofIm = DofMap.getCompDoF(mortar.nodesMaster);
dofM = DofMap.getCompDoF(1:masterMesh.nNodes);
dofM = dofM(~ismember(dofM,dofIm));
dofIs = DofMap.getCompDoF(mortar.nodesSlave);
dofS = DofMap.getCompDoF(1:slaveMesh.nNodes);
dofS = dofS(~ismember(dofS,dofIs));
dofMult = dofIs;
Kmm = KMaster(dofM,dofM);
KmIm = KMaster(dofM,dofIm);
Kss = KSlave(dofS,dofS);
KsIs = KSlave(dofS, dofIs);
KImIm = KMaster(dofIm,dofIm);
KIsIs = KSlave(dofIs,dofIs);
% expanding the mortar operator also to the y direction
Erbf = Drbf\Mrbf;
Esb = Dsb\Msb;
switch integration
   case 'RBF'
      D = expandMat(Drbf,2);
      M = expandMat(Mrbf,2);
   case 'SB'
      D = expandMat(Dsb,2);
      M = expandMat(Msb,2);
   case 'P0'
      D = expandMat(Dconsistent,2);
      M = expandMat(Mconst,2);
end

% removing endpoint nodal multipliers
if bound && ~strcmp(integration,'P0')
   D([5 6],:) = D([5 6],:) + D([1 2],:);
   D([end end-1],:) = D([end end-1],:) + D([3 4],:);
   M([5 6],:) = M([5 6],:) + M([1 2],:);
   M([end end-1],:) = M([end end-1],:) + M([3 4],:);
   D([1 2 3 4],:) = [];
   M([1 2 3 4],:) = [];
   dofMult = DofMap.getCompDoF(mortar.nodesSlave(3:end),2);
   nMult = numel(dofMult);
   H = zeros(nMult,nMult);   
else
   dofMult = DofMap.getCompDoF(mortar.nodesSlave(1:end),2);
end

switch integration 
   case 'P0'
      dofMult = DofMap.getCompDoF((1:mortar.nElSlave)');
      nMult = numel(dofMult);
      H1 = computeStabilizationMatrix(mortar,Em,Es);
      [Dp0,Mp0,Dcp0] = mortar.computeMortarConstant(nG,4);
      H3 = computeStabilizationMatrix3(mortar,Dp0,Mp0,KMaster,KSlave);
      %H4 = computeStabilizationMatrix4(mortar,Dp0,Mp0,KMaster,KSlave);
      H = H3;
end
% complete saddle point matrix
r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofMult))];
r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofMult))] ;
r3 = [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M'];
r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D'];
r5 = [zeros(length(dofMult), length(dofM)), zeros(length(dofMult), length(dofS)), -M, D, -H];
K = [r1;r2;r3;r4;r5];


% Fix displacements BCS - normal component of boundary nodes
if bound || strcmp(integration,'P0')
   dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
   dirDoFSlaveX = dof.getDoF(dirNod,'slave',2);
   dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 3,:));
   dirDoFSlaveY = dof.getDoF(dirNod,'slave',2);
   %
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
   dirDoFMasterX = dof.getDoF(dirNod,'master',2);
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
   dirDoFMasterY = dof.getDoF(dirNod,'master',2);
   dirDoF = unique([dirDoFMasterX;dirDoFMasterY;dirDoFSlaveX;dirDoFSlaveY]);
else
   dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
   dirNod = dirNod(~ismember(dirNod,mortar.nodesSlave));
   dirDoFSlaveX = dof.getDoF(dirNod,'slave',2);
   dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 3,:));
   dirDoFSlaveY = dof.getDoF(dirNod,'slave',2);
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
   dirNod = dirNod(~ismember(dirNod,mortar.nodesMaster));
   dirDoFMasterX = dof.getDoF(dirNod,'master',2);
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
   dirDoFMasterY = dof.getDoF(dirNod,'master',2);
   dirDoF = unique([dirDoFMasterX;dirDoFMasterY;dirDoFSlaveX;dirDoFSlaveY]);
end
% K = full(K);
K(dirDoF,:) = [];
K(:,dirDoF) = [];

% Inf-sup test: Elman 
% tic 
A = K(1:end-nMult,1:end-nMult);
B = K(end-nMult+1:end,1:end-nMult);
hS = 1/NSX;
X = A\B';
S = B*X;
S = full(S);
q = (sum(D,2));
%H = zeros(nMult,nMult);  
invQ = full(diag(1./(hS*q))); % scaled by the length for proper dual norm bound
eS = eig(invQ*(S+H));
infSup(iref) = sqrt(min(eig(invQ*(S+H))));
end

%% inf-sup constant for stabilized formulation
% Plotting the results of the inf-sup test on the macroelement
h = [1/4; 1/8; 1/16; 1/32; 1/64];
out = load("outStab.mat");

out1 = out.case1;
out2 = out.case2;
out3 = out.case3;

% Create figure
figure(1)
hold on
grid on

% Plot coarse mesh data
plot(log(h), out1, 'k-s', 'LineWidth', 1, 'DisplayName', 'r_h = 1/2$');
hold on 
plot(log(h), out2, 'k-s', 'LineWidth', 1, 'DisplayName', 'r_h = 1/3$');
plot(log(h), out3, 'k-s', 'LineWidth', 1, 'DisplayName', 'r_h = 1/3$');

% Labels and title with LaTeX
xlabel('$\log(h)$', 'Interpreter', 'latex');
ylabel('$\beta_h^{\star}$', 'Interpreter', 'latex');
% legend.show()
%title('Inf-Sup Test Results', 'Interpreter', 'latex');
%ylim([-0.3 0.1])

xlabel('$\log(h)$', 'Interpreter', 'latex');
ylabel('$\log(\beta_h^{\star})$', 'Interpreter', 'latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;


% Export high-resolution figure
exportgraphics(figure(1), 'infSup_stabilized.pdf', 'Resolution', 300);

