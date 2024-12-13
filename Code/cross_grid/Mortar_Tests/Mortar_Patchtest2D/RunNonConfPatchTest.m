function [mult_x,mult_y] = RunNonConfPatchTest(mult_type,Fx,Fy,Dmat)
% Run a 2 non conforming block patch test with standard or dual lagrange
% multipliers

% plots solution to paraview
% output: multipliers on the interface

% DEFINE MODEL

% IMPORT MESHES
topMesh = Mesh();
bottomMesh = Mesh();


% gauss class for 1D element integration
nG = 6;

fileNameTop = 'Mesh_flat/TopBlock_hexa.msh';
fileNameBottom = 'Mesh_flat/BottomBlock_hexa.msh';

% Set the input file name
gaussQuad = Gauss(12,2,2);

flagTop = 'slave';

integration = "RBF"; % RBF or SB
solution_scheme = "SP"; % SP (saddle point) % COND (condensated)


% Import the mesh data into the Mesh object
topMesh.importGMSHmesh(fileNameTop);
bottomMesh.importGMSHmesh(fileNameBottom);

% assign domains to mortar or slave tag
% decide if the top domain is master or slave

if strcmp(flagTop, 'master')
   bottom = 'slave';
   masterMesh = topMesh;
   slaveMesh = bottomMesh;
elseif strcmp(flagTop, 'slave')
   bottom = 'master';
   slaveMesh = topMesh;
   masterMesh = bottomMesh;
end

% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh,gaussQuad);
elemsSlave = Elements(slaveMesh,gaussQuad);

% computing stiffness matrix on each domain
% DOMAIN 1
KMaster = stiff(masterMesh, elemsMaster, Dmat, gaussQuad);
KSlave = stiff(slaveMesh, elemsSlave, Dmat, gaussQuad);
% get id of nodes belonging to master and slave interfaces
nodesMaster = unique(masterMesh.edges(masterMesh.edgeTag == 1,:));
nodesSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 1,:));

mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
% compute mortar operator
if strcmp(integration,'RBF')
   [Dtmp, Mtmp] = mortar.computeMortarRBF(nG,4,'gauss',mult_type);
elseif strcmp(integration, 'EB')
   [Dtmp, Mtmp] = mortar.computeMortarElementBased(nG);
elseif strcmp(integration, 'SB')
   [Dtmp, Mtmp] = mortar.computeMortarSegmentBased(3,mult_type);
end



Etmp = Dtmp\Mtmp;
% reordering the matrix of the system
%
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
% get Loaded dofs on top edge
nodesLoad = unique(topMesh.edges(topMesh.edgeTag == 2,:));
% special treatment of infextreme points (having force 1/2)
n_ext = nodesLoad(ismember(topMesh.coordinates(nodesLoad,1),[0; 1]));
loadDoFext = 2*n_ext; % loading y component
loadDoFY = 2*nodesLoad(~ismember(nodesLoad, n_ext)); % loading y component
loadDoFext = getGlobalDofs(loadDoFext, dofM,dofS,dofIm,dofIs, flagTop);
loadDoFY = getGlobalDofs(loadDoFY, dofM,dofS,dofIm,dofIs, flagTop);
f(loadDoFext) = Fy*(0.5/(length(nodesLoad)-1));
f(loadDoFY) = Fy*(1/(length(nodesLoad)-1));

%-------------------------LATERAL LOAD BCS ------------------------
nodesLoad = find(abs(topMesh.coordinates(:,1)-0)<1e-3);
% get top node (it has half of the entities influence)
[yC,id] = sort(topMesh.coordinates(nodesLoad,2),'ascend');
l1 = abs(yC(end)-yC(end-1));
l2 = abs(yC(2)-yC(1));
n1 = nodesLoad(id(end));
f(getGlobalDofs(2*n1-1,dofM,dofS,dofIm,dofIs,flagTop)) = Fx*l1/2;
nIn = nodesLoad(id(3:end-1));
f(getGlobalDofs(2*nIn-1,dofM,dofS,dofIm,dofIs,flagTop)) = Fx*l1;
n2 = nodesLoad(id(2));
f(getGlobalDofs(2*n2-1,dofM,dofS,dofIm,dofIs,'slave')) = Fx*(l1/2+l2/2);
n3 = nodesLoad(id(1));
f(getGlobalDofs(2*n3-1,dofM,dofS,dofIm,dofIs,'interfaceSlave')) = Fx*l2/2;
% 
% % LATERAL LOAD ALSO TO MASTER SIDE (BOTTOM)
% % get top node (it has half of the entities influence)
% nodesLoad = find(abs(bottomMesh.coordinates(:,1)-0)<1e-3);
% [yC,id] = sort(bottomMesh.coordinates(nodesLoad,2),'descend');
% l = abs(yC(end)-yC(end-1));
% n1 = nodesLoad(id(end));
% f(getGlobalDofs(2*n1-1,dofM,dofS,dofIm,dofIs,'interfaceMaster')) = Fx*l/2;
% nIn = nodesLoad(id(2:end-1));
% f(getGlobalDofs(2*nIn-1,dofM,dofS,dofIm,dofIs,'master')) = Fx*l;
% n2 = nodesLoad(id(1));
% f(getGlobalDofs(2*n2-1,dofM,dofS,dofIm,dofIs,'master')) = Fx*l/2;

%------------------- BOTTOM FIXED BCS -----------------------------
% get fixed dofs on bottom edge
% y bottom constraint
dirNod = find(abs(bottomMesh.coordinates(:,2)-0)<1e-3);
dirBotDoF = getGlobalDofs(2*dirNod,dofM,dofS,dofIm,dofIs,bottom);
[K,f(1:nf)] = applyDir(dirBotDoF, zeros(length(dirBotDoF),1), K, f(1:nf));

% -------------------LATERAL CONSTRAINT BCS-----------------------
% get nodes on right edge of master domain
nodesFixX = find(abs(topMesh.coordinates(:,1)-1)<1e-3);
dirXDoF = getGlobalDofs(2*nodesFixX-1,dofM,dofS,dofIm,dofIs,'slave');
id = abs(topMesh.coordinates(nodesFixX,2)-1)<1e-3;
dirXIntDoF = getGlobalDofs(2*nodesFixX(id)-1,dofM,dofS,dofIm,dofIs,'interfaceSlave');
[K,f(1:nf)] = applyDir(dirXDoF, zeros(length(dirXDoF),1), K, f(1:nf));
%[K,f(1:nf)] = applyDir(dirXIntDoF, zeros(length(dirXIntDoF),1), K, f(1:nf));

% get nodes on right edge of master domain
nodesFixX = find(abs(bottomMesh.coordinates(:,1)-1)<1e-3);
dirXDoF = getGlobalDofs(2*nodesFixX-1,dofM,dofS,dofIm,dofIs,'master');
id = abs(bottomMesh.coordinates(nodesFixX,2)-1)<1e-3;
dirXIntDoF = getGlobalDofs(2*nodesFixX(id)-1,dofM,dofS,dofIm,dofIs,'interfaceMaster');
[K,f(1:nf)] = applyDir(dirXDoF, zeros(length(dirXDoF),1), K, f(1:nf));
[K,f(1:nf)] = applyDir(dirXIntDoF, zeros(length(dirXIntDoF),1), K, f(1:nf));

n = numel(dofM)+numel(dofS);

% modify rhs array for condensated system
if strcmp(solution_scheme,"COND")
   f(n+1:n+numel(dofIm)) = ...
      f(n+1:n+numel(dofIm)) + ...
      E'*f(n+numel(dofIm)+1:n+numel(dofIm)+numel(dofIs));
end

fIntSlave = f(n+numel(dofIm)+1:n+numel(dofIm)+numel(dofIs));
if strcmp(solution_scheme,"COND")
   f = f(1:n+numel(dofIm)); % remove entries
end

% solve linear system
u = K\f;

n = length(dofM)+length(dofS);
l = length(dofM)+length(dofS)+length(dofIm);
% get slave side solution
if strcmp(solution_scheme, "SP")
   u_slave = u(l+1:l+length(dofIs));
   mult = u(l+numel(dofIs)+1:end);
else
   u_slave = E*u(length(dofM)+length(dofS)+1:end);
   mult = (D')\(fIntSlave-KsIs'*u(numel(dofM)+1:n)-KIsIs*E*u(n+1:l)); % recover multipliers
end
% plot solution (use the plot_function already used for the RBF stand alone
% tests!)

u_top = zeros(topMesh.nNodes,1);
u_bottom = zeros(bottomMesh.nNodes,1);

%collect displacement of master domain and slave domain, according to user assignment;
if strcmp(flagTop, 'master')
   u_top(dofM) = u(1:length(dofM));
   u_top(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
   u_bottom(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
   u_bottom(dofIs) = u_slave;
elseif strcmp(flagTop, 'slave')
   u_bottom(dofM) = u(1:length(dofM));
   u_bottom(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
   u_top(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
   u_top(dofIs) = u_slave;
end

% compute stresses
stressMaster = computeStress(masterMesh,elemsMaster,Dmat,u_bottom,gaussQuad);
stressSlave = computeStress(slaveMesh,elemsSlave,Dmat,u_top,gaussQuad);
%
plotSolution(topMesh,strcat(mult_type,'_top'),u_top,stressSlave);
plotSolution(bottomMesh,strcat(mult_type,'_bottom'),u_bottom,stressMaster);
%
% get stress on slave interface
cellInterf = any(ismember(slaveMesh.surfaces,nodesSlave),2);
s_x = stressSlave(cellInterf,3);
s_y = stressSlave(cellInterf,2);

% plot lagrange multipliers in matlab
% compute multipliers
[~,id] = sort(slaveMesh.coordinates(mortar.nodesSlave,1));
mult_x = mult(2*id-1);
mult_y = mult(2*id);
end

