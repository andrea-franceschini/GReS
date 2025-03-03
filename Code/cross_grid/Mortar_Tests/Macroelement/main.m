clear; close all; 

%% MacroElement approach to spot violation of inf-sup condition

% Input: specify number of elements for top (master) and bottom (slave)
% side of the mesh.

% minimum number of elements is 2 for both sides

% all boundary nodes have fixed displacements.
% each element has fixed size

nG = 6;
E = 1; 
nu = 0.25;
integration = "P0"; % RBF,SB,P0

Dmat = zeros(3);
Dmat([1 5]) = 1-nu;
Dmat([2 4]) = nu;
Dmat(9) = 0.5*(1-2*nu);
Dmat = (E/((1+nu)*(1-2*nu)))*Dmat;

% Set the input file name
gaussQuad = Gauss(12,2,2);

NM0 = 2;
rat = 2;
nR = 1;
h = zeros(nR,1);
inf = zeros(nR,1);
sup = inf;


for iref = 1:nR
% read mesh
NM = NM0*2^(iref-1);
NS = rat*NM;
h(iref) = 1/NS;
getPatchMesh('Mesh/master.geo','master',NM);
getPatchMesh('Mesh/slave.geo','slave',NS);
[slaveMesh,masterMesh] = deal(Mesh(),Mesh()); 
slaveMesh.importGMSHmesh('Mesh/slave.msh')
masterMesh.importGMSHmesh('Mesh/master.msh');

% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh,gaussQuad);
elemsSlave = Elements(slaveMesh,gaussQuad);


KMaster = stiff(masterMesh, elemsMaster, Dmat, gaussQuad);
KSlave = stiff(slaveMesh, elemsSlave, Dmat, gaussQuad);


mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
dof = DofMap(masterMesh,mortar.nodesMaster,slaveMesh,mortar.nodesSlave);


% compute mortar operator1
[Drbf, Mrbf] = mortar.computeMortarRBF(nG,4,'gauss','standard');
[Dsb, Msb] = mortar.computeMortarSegmentBased(nG,'standard');
[Dconst,Mconst,Dconsistent] = mortar.computeMortarConstant(nG,4);

dofIm = DofMap.getCompDoF(mortar.nodesMaster);
dofM = DofMap.getCompDoF(1:masterMesh.nNodes);
dofM = dofM(~ismember(dofM,dofIm));
dofIs = DofMap.getCompDoF(mortar.nodesSlave);
dofS = DofMap.getCompDoF(1:slaveMesh.nNodes);
dofS = dofS(~ismember(dofS,dofIs));
dofMult = dofIs;
switch integration 
   case 'P0'
      dofMult = DofMap.getCompDoF((1:mortar.nElSlave)');
      nMult = numel(dofMult);
end
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
      H = mortar.computePressureJumpMat();
      H = (1/NS)*expandMat(H,2);
end

if ~strcmp(integration,'P0')
   D([5 6],:) = D([5 6],:) + D([1 2],:);
   D([end end-1],:) = D([end end-1],:) + D([3 4],:);
   M([5 6],:) = M([5 6],:) + M([1 2],:);
   M([end end-1],:) = M([end end-1],:) + M([3 4],:);
   D([1 2 3 4],:) = [];
   M([1 2 3 4],:) = [];
   dofMult = DofMap.getCompDoF(mortar.nodesSlave(3:end),2);
   nMult = numel(dofMult);
   H = zeros(nMult,nMult);
end

% complete saddle point matrix
r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofMult))];
r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofMult))] ;
r3 = [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M'];
r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D'];
r5 = [zeros(length(dofMult), length(dofM)), zeros(length(dofMult), length(dofS)), -M, D, -H];
K = [r1;r2;r3;r4;r5];


% Fix displacements BCS
dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
dirDoFSlave = dof.getDoF(dirNod,'slave',2);
dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
dirDoFMaster = dof.getDoF(dirNod,'master',2);
dirDoF = [dirDoFMaster;dirDoFSlave];
% K = full(K);
K(dirDoF,:) = [];
K(:,dirDoF) = [];

%u = K\f;

% Eigen-decomposition of Schur complement matrix
A = K(1:end-nMult,1:end-nMult);
B = K(end-nMult+1:end,1:end-nMult);
tic 
X = A\B';
S = B*X;

% extract only some row/cols from B
id = sum(full(B) ~= 0,1)==0;
Bred = B(:,~id);
Ared = A(~id,~id);
Xred = Ared\Bred';
Sred = Bred*Xred;

% Bt = B;
% Bt(:,[1 2]) = B(:,[3 4]); 
% Bt(:,[3 4]) = B(:,[1 2]);
% Bt = -Bt;
% H1 = Bt'*Bt;

% B(:,)
[vS,eS] = eig(full(S));
inf(iref) = sqrt(min(abs(diag(eS))));
sup(iref) = sqrt(max(abs(diag(eS))));
end
%S = 0.5*(S+S');
loglog(h, inf, 'k-o') % Grafico log-log
hold on

for i = 1:length(inf)-1
    ratio = inf(i+1) / inf(i); % Calcola il rapporto tra due valori consecutivi
    text(h(i+1), inf(i+1), sprintf('%.2f', ratio), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% [vs,lambdaS] = eig(S);
% [vh,lambdaH] = eig(full(H));
% e1 = -((1+nu)*(1-2*nu)/(E*(1-nu)));
% e2 = 3/(2+(1-2*nu)/(1-nu));
% e = e1*e2;
% fprintf('Condition number of Schur complement matrix: %2.5e \n',cond(S));
% fprintf('Smallest eigenvalue: %2.5e \n',min(diag(abs(lambdaS))));