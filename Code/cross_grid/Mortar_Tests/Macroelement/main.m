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

NM = 2;
NS = 2;
% read mesh
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


% compute mortar operator1.2
[Drbf, Mrbf] = mortar.computeMortarRBF(nG,4,'gauss','dual');
[Dsb, Msb] = mortar.computeMortarSegmentBased(nG,'dual');
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
end
% complete saddle point matrix
r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofMult))];
r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofMult))] ;
r3 = [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M'];
r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D'];
r5 = [zeros(length(dofMult), length(dofM)), zeros(length(dofMult), length(dofS)), -M, D, zeros(numel(dofMult),numel(dofMult))];
K = [r1;r2;r3;r4;r5];


% Fix displacements BCS
dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
dirDoFSlave = dof.getDoF(dirNod,'slave',2);
dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
dirDoFMaster = dof.getDoF(dirNod,'master',2);
dirDoF = [dirDoFMaster;dirDoFSlave];
K = full(K);
K(dirDoF,:) = [];
K(:,dirDoF) = [];

nMult = numel(dofMult);
% eliminate multipliers on fault tip (only for nodal multipliers)
if ~strcmp(integration,'P0')
   % sum columns in B^t
   K(:,end-numel(dofMult)+5) = K(:,end-numel(dofMult)+5)+K(:,end-numel(dofMult)+1);
   K(:,end-1) = K(:,end-1)+K(:,end-numel(dofMult)+3);
   K(:,end-numel(dofMult)+6) = K(:,end-numel(dofMult)+6)+K(:,end-numel(dofMult)+2);
   K(:,end) = K(:,end)+K(:,end-numel(dofMult)+4);
   %sum rows in B
   K(end-numel(dofMult)+5,:) = K(end-numel(dofMult)+5,:)+K(end-numel(dofMult)+1,:);
   K(end-1,:) = K(end-1,:)+K(end-numel(dofMult)+3,:);
   K(end-numel(dofMult)+6,:) = K(end-numel(dofMult)+6,:)+K(end-numel(dofMult)+2,:);
   K(end,:) = K(end,:)+K(end-numel(dofMult)+4,:);
   % delete multipliers on the tip
   K(end-numel(dofMult)+1,:) = [];
   K(end-numel(dofMult)+2,:) = [];
   K(end-numel(dofMult)+3,:) = [];
   K(end-numel(dofMult)+4,:) = [];
   K(:,end-numel(dofMult)+1) = [];
   K(:,end-numel(dofMult)+2) = [];
   K(:,end-numel(dofMult)+3) = [];
   K(:,end-numel(dofMult)+4) = [];
   nMult = nMult-4;
end


% Eigen-decomposition of Schur complement matrix
A = K(1:end-nMult,1:end-nMult);
B = K(end-nMult+1:end,1:end-nMult);
S = -B*inv(A)*B';
S = 0.5*(S+S');
[v,lambda] = eig(S);

e1 = -((1+nu)*(1-2*nu)/(E*(1-nu)));
e2 = 3/(2+(1-2*nu)/(1-nu));
e = e1*e2;
fprintf('Condition number of Schur complement matrix: %2.5e \n',cond(S));
fprintf('Smallest eigenvalue: %2.5e \n',min(diag(abs(lambda))));