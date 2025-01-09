function [mult_x,mult_y,u_slave,u_master] = RunNonConfPatchTest(mult_type,Fx,Fy,DmatS,DmatM,patch,stab,nXs,nYs,nXm,nYm)
% Run a 2 non conforming block patch test with standard or dual lagrange
% multipliers

% plots solution to paraview
% output: multipliers on the interface

% gauss class for 1D element integration
nG = 16;

% Set the input file name
gaussQuad = Gauss(12,2,2);

integration = "SB"; % RBF or SB
solution_scheme = "SP"; % SP (saddle point) % COND (condensated)

masterMesh = getMesh('Mesh/bottomBlock.geo','bottom',nXm,nYm);
slaveMesh = getMesh('Mesh/topBlock.geo','top',nXs,nYs);
% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh,gaussQuad);
elemsSlave = Elements(slaveMesh,gaussQuad);

% computing stiffness matrix on each domain
% DOMAIN 1
KMaster = stiff(masterMesh, elemsMaster, DmatM, gaussQuad);
KSlave = stiff(slaveMesh, elemsSlave, DmatS, gaussQuad);
% get id of nodes belonging to master and slave interfaces
nodesMaster = unique(masterMesh.edges(masterMesh.edgeTag == 1,:));
nodesSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 1,:));


dof = DofMap(masterMesh,nodesMaster,slaveMesh,nodesSlave);
% % remove extreme nodes in contact with the boundary
mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
% compute mortar operator
if strcmp(integration,'RBF')
   [Dtmp, Mtmp] = mortar.computeMortarRBF(2,4,'gauss',mult_type);
elseif strcmp(integration, 'EB')
   [Dtmp, Mtmp] = mortar.computeMortarElementBased(nG);
elseif strcmp(integration, 'SB')
   [Dtmp, Mtmp] = mortar.computeMortarSegmentBased(2,mult_type);
end

if(strcmp(mult_type,'dual'))
   Dtmp = diag(sum(Dtmp,2));
end


Etmp = Dtmp\Mtmp;
% reordering the matrix of the system
%
%Etmp(abs(Etmp)<1e-8) = 0;
% | dofM = nodi interni master      |
% | dofS = nodi interni slave       |
% | dofIm = nodi interfaccia master |
% | dofIs = nodi interfaccia slave  |
% dof Lagrange multiplier not needed

% build an easy to access dof map
dofIm = getDoF(nodesMaster');
dofM = getDoF(1:masterMesh.nNodes);
dofM = dofM(~ismember(dofM,dofIm));
dofIs = getDoF(nodesSlave');
dofS = getDoF(1:slaveMesh.nNodes);
dofS = dofS(~ismember(dofS,dofIs));
Kmm = KMaster(dofM,dofM);
KmIm = KMaster(dofM,dofIm);
Kss = KSlave(dofS,dofS);
KsIs = KSlave(dofS, dofIs);
KImIm = KMaster(dofIm,dofIm);
KIsIs = KSlave(dofIs,dofIs);
% expanding the mortar operator also to the y direction
nM = length(nodesMaster);
nS = length(nodesSlave);
E = zeros(2*nS,2*nM);
M = zeros(2*nS,2*nM);
D = zeros(2*nS,2*nS);
E(2*(1:nS)'-1,2*(1:nM)'-1) = Etmp;
E(2*(1:nS)',2*(1:nM)') = Etmp;
M(2*(1:nS)'-1,2*(1:nM)'-1) = Mtmp;
M(2*(1:nS)',2*(1:nM)') = Mtmp;
D(2*(1:nS)'-1,2*(1:nS)'-1) = Dtmp;
D(2*(1:nS)',2*(1:nS)') = Dtmp;

% condensated stiffness matrix
KCOND = [Kmm, zeros(length(dofM),length(dofS)), KmIm;
   zeros(length(dofS),length(dofM)), Kss, KsIs*E;
   KmIm', E'*KsIs', KImIm+E'*KIsIs*E];

% complete saddle point matrix
KSP = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofIs));
   zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofIs)) ;
   KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M';
   zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D';
   zeros(length(dofIs), length(dofM)), zeros(length(dofIs), length(dofS)), -M, D,  zeros(length(dofIs), length(dofIs))
   ];

% creating global dof array and initializing forcing vector
listDofsCOND = [dofM';dofS';dofIm'];
listDofsSP = [dofM';dofS';dofIm'; dofIs'; dofIs'];
nCond = numel(listDofsCOND);
nSP = numel(listDofsSP);
f = zeros(length(listDofsSP),1);

if strcmp(solution_scheme, "SP")
   K = KSP;
   nf = nSP;
else
   K = KCOND;
   nf = nCond;
end
% ------------------- APPLY BCS -------------------------------

%------------------- TOP LOAD BCS -----------------------------

if patch == 1
   % get Loaded dofs on top edge
   [nodesLoad,lInf] = getAreaInf(slaveMesh,2);
   % special treatment of extreme points (having force 1/2)
   loadDoF = dof.getDoF(nodesLoad,'slave','y');
   f(loadDoF) = Fy*lInf;
end

%-------------------------LATERAL LOAD BCS ------------------------
if patch == 2 || patch == 3
   [nodesLoad,lInf] = getAreaInf(slaveMesh,3);
   loadDoF = dof.getDoF(nodesLoad,'slave','x');
   f(loadDoF) = Fx*lInf;
end


%------------------- BOTTOM FIXED BCS -----------------------------
% get fixed dofs on bottom edge
% y bottom constraint
dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
if patch == 1 || patch == 2
   dirBotDoF = dof.getDoF(dirNod,'master');
elseif patch == 3
   dirBotDoF = dof.getDoF(dirNod,'master','y');
end
[K,f(1:nf)] = applyDir(dirBotDoF, zeros(length(dirBotDoF),1), K, f(1:nf));

%------------------- LATERAL FIXED BCS -----------------------------
if patch == 3
   dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 4,:));
   % remove interface slave node for unstabilized version
   if strcmp(stab,'unstable')
      [~,id] = min(slaveMesh.coordinates(dirNod,2));
      dirNod(id) = [];
   end
   dirDoF = dof.getDoF(dirNod,'slave','x');
   [K,f(1:nf)] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f(1:nf));
   %
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
   dirDoF = dof.getDoF(dirNod,'master','x');
   [K,f(1:nf)] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f(1:nf));
end

masterDof = dof.getDoF(nodesMaster(3:end),'master');
slaveDof = dof.getDoF(nodesSlave(3:end),'slave');
lagDof = dof.getDoF(nodesSlave(3:end),'lagrange');
if numel(slaveDof)<20
A = K([masterDof;slaveDof],[masterDof;slaveDof]);
B = K([masterDof;slaveDof],lagDof);
S = -(B')*(inv(A)*B);
[V,Eigen] = eig(full(S));
end
% remove lagrange multipliers belonging to the interface (then use constant
% interpolation and extend the adjacent value)
% get boundary nodes and adjacent node id
if strcmp(stab,'stable')
   [K,f,dofs] = handleEndPoints(slaveMesh,nodesSlave,dof,K,f);
end

% solve linear system
u = K\f;

if strcmp(stab,'stable')
   u = addSolToEndPoints(slaveMesh,nodesSlave,u,dof);
end

n = length(dofM)+length(dofS);
l = length(dofM)+length(dofS)+length(dofIm);
% get slave side solution
u_master = u(n+1:n+length(dofIm));
u_slave = u(l+1:l+length(dofIs));
mult = u(l+numel(dofIs)+1:end);
% plot solution (use the plot_function already used for the RBF stand alone
% tests!)

u_top = zeros(slaveMesh.nNodes,1);
u_bottom = zeros(masterMesh.nNodes,1);

%collect displacement of master domain and slave domain, according to user assignment;
u_bottom(dofM) = u(1:length(dofM));
u_bottom(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
u_top(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
u_top(dofIs) = u_slave;

% compute stresses
stressMaster = computeStress(masterMesh,elemsMaster,DmatM,u_bottom,gaussQuad);
stressSlave = computeStress(slaveMesh,elemsSlave,DmatS,u_top,gaussQuad);
%
plotSolution(slaveMesh,strcat(mult_type,'_top'),u_top,stressSlave);
plotSolution(masterMesh,strcat(mult_type,'_bottom'),u_bottom,stressMaster);
%
% get stress on slave interface
% cellInterf = any(ismember(slaveMesh.surfaces,nodesSlave),2);
% s_x = stressSlave(cellInterf,3);
% s_y = stressSlave(cellInterf,2);

% lagrange multipliers post processing
% [Htmp,h] = mortar.computeWeighMat(mult_type);
% nh = numel(h);
% H = zeros(2*nh,2*nh);
% h = repelem(h,2,1);
% H(2*(1:nh)'-1,2*(1:nh)'-1) = Htmp;
% H(2*(1:nh)',2*(1:nh)') = Htmp;
% H = H./h;
%mult = H*mult;

% invert master and slave for post-processing
mortarSwap = Mortar2D(1,slaveMesh,1,masterMesh,1);
if strcmp(integration,'RBF')
   [Dswap, Mtmp] = mortarSwap.computeMortarRBF(2,4,'gauss',mult_type);
elseif strcmp(integration, 'EB')
   [Dswap, Mtmp] = mortarSwap.computeMortarElementBased(nG);
elseif strcmp(integration, 'SB')
   [Dswap, Mtmp] = mortarSwap.computeMortarSegmentBased(2,mult_type);
end

if(strcmp(mult_type,'dual'))
   Dswap = diag(sum(Dswap,2));
end
Eswap = Dswap\Mtmp;


% plot lagrange multipliers in matlab
% compute multipliers
[~,id] = sort(masterMesh.coordinates(mortar.nodesMaster,1));
u_master = [u_master(2*id-1) u_master(2*id)];


[~,id] = sort(slaveMesh.coordinates(mortar.nodesSlave,1));
mult_x = mult(2*id-1);
mult_y = mult(2*id);
u_slave = [u_slave(2*id-1) u_slave(2*id)];
Esw = zeros(2*nM,2*nS);
Esw(2*(1:nM)'-1,2*(1:nS)'-1) = Eswap;
Esw(2*(1:nM)',2*(1:nS)') = Eswap;
%mult = E*(Esw*mult);
mult_x = mult(2*id-1);
mult_y = mult(2*id);

end

