%clear; close all; 

%% MacroElement approach to spot violation of inf-sup condition

% Input: specify number of elements for top (master) and bottom (slave)
% side of the mesh.

% minimum number of elements is 2 for both sides

% all boundary nodes have fixed displacements.
% each element has fixed size

nG = 6;
E = 10000; 
nu = 0.25;
integration = "P0"; % RBF,SB,P0
mult_type = 'standard';
bound = false;


Dmat = zeros(3);
Dmat([1 5]) = 1-nu;
Dmat([2 4]) = nu;
Dmat(9) = 0.5*(1-2*nu);
Dmat = (E/((1+nu)*(1-2*nu)))*Dmat;

% Set the input file name
gaussQuad = Gauss(12,2,2);

NM0X = 5;
NS0X = 5;
NM0Y = 1;
NS0Y = 1;
nR = 1;
h = zeros(nR,1);
beta1 = zeros(nR,1);
beta2 = zeros(nR,1);
beta3 = zeros(nR,1);

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

% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh,gaussQuad);
elemsSlave = Elements(slaveMesh,gaussQuad);


KMaster = stiff(masterMesh, elemsMaster, Dmat, gaussQuad);
KSlave = stiff(slaveMesh, elemsSlave, Dmat, gaussQuad);


mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
dof = DofMap(masterMesh,mortar.nodesMaster,slaveMesh,mortar.nodesSlave);

switch integration
   case 'SB'
      P = mortar.computeMassMat(mult_type);
   case 'P0'
      P = mortar.computeMassMatP0();
end
P = expandMat(P,2);

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
      H = mortar.computePressureJumpMat();
      H = (1/NSX)*expandMat(H,2);
end

% removing endpoint nodal multipliers
if bound && ~strcmp(integration,'P0')
   D([5 6],:) = D([5 6],:) + D([1 2],:);
   D([end end-1],:) = D([end end-1],:) + D([3 4],:);
   M([5 6],:) = M([5 6],:) + M([1 2],:);
   M([end end-1],:) = M([end end-1],:) + M([3 4],:);
   D([1 2 3 4],:) = [];
   M([1 2 3 4],:) = [];
   P([1 2 3 4],:) = [];
   P(:,[1 2 3 4]) = [];
   dofMult = DofMap.getCompDoF(mortar.nodesSlave(3:end),2);
else
   dofMult = DofMap.getCompDoF(mortar.nodesSlave(1:end),2);
end
nMult = numel(dofMult);
switch integration 
   case 'P0'
      dofMult = DofMap.getCompDoF((1:mortar.nElSlave)');
      nMult = numel(dofMult);
end
H = zeros(nMult,nMult);
%H = zeros(nMult,nMult);
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

%u = K\f;

% Inf-sup test : BATHE-CHAPELLE
A = K(1:end-nMult,1:end-nMult);
B = K(end-nMult+1:end,1:end-nMult);
hS = 1/NSX;

PB = P\B;
G = (1/hS)*(B'*PB);
% solve the generalized eigenvalue problem computing only necessary eigs
e_sparse = real(eigs(G,A,size(B,1),'la'));

id = abs(e_sparse)<1e-13;
k = sum(id)+size(G,1)-size(B,1);
m = size(B,2); n = size(B,1);
kpm = k-(m-n); % compute the number of spurious pressure modes
e = sqrt(min(e_sparse(~id)));
beta1(iref) = e;

% Inf-sup test: trying using proper fractional norms
% A = K(1:end-nMult,1:end-nMult);
% B = K(end-nMult+1:end,1:end-nMult);
% tic 
X = A\B';
S = B*X;
% eS = min(real(eig(full(S))));
% eA = max(real(eig(full(A))));
% beta3(iref) = sqrt(eS/eA);



% % Inf-sup test : BOFFI-BREZZI-FIRTIN
% svd of matrix B modified with suitable norm matrices...
% smallest positive singular value of Sy*B*Sy where
% [V, D] = eigs(P, size(P,1));
% Sy = (hS^(-1/4))*(V * sqrt(D) * V');
% [V, D] = eigs(A, size(A,1));
% 
% Sx = V * sqrt(D) * V';
% e2 = svd(Sy*(B*Sx));
% beta2(iref) = min(e2);
% Sy = sqrtm(P);
% Sx = sqrt(A);
plot(log(h),log(beta1), 'k-o') 
hold on
% plot(log(h),log(beta2), 'r-o')
% plot(log(h),log(beta3), 'b-o')
legend('Bathe-Chapelle','Brezzi')
end

%%
% for i = 1:length(inf)-1
%     ratio = inf(i+1) / inf(i); % Calcola il rapporto tra due valori consecutivi
%     text(h(i+1), inf(i+1), sprintf('%.2f', ratio), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end

% [vs,lambdaS] = eig(S);
% [vh,lambdaH] = eig(full(H));
% e1 = -((1+nu)*(1-2*nu)/(E*(1-nu)));
% e2 = 3/(2+(1-2*nu)/(1-nu));
% e = e1*e2;
% fprintf('Condition number of Schur complement matrix: %2.5e \n',cond(S));
% fprintf('Smallest eigenvalue: %2.5e \n',min(diag(abs(lambdaS))));