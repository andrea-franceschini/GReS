function [mult_x,mult_y,u_slave,u_master] = RunConfPatch(Fx,Fy,DmatS,DmatM,patch,stab,nX,nY,nYm)
% Run a 2 non conforming block patch test with standard or dual lagrange
% multipliers

% plots solution to paraview
% output: multipliers on the interface

% DEFINE MODEL

% IMPORT MESHES
masterMesh = getMesh('Mesh/bottomBlock.geo','bottomConf',nX,nYm);
slaveMesh = getMesh('Mesh/topBlock.geo','topConf',nX,nY);

% Set the input file name
gaussQuad = Gauss(12,2,2);


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
% compute cross grid conforming operator
[HMtmp,HDtmp] = mortar.computeConfCrossGridMat();

% reordering the matrix of the system
%
% | dofM = nodi interni master      |
% | dofS = nodi interni slave       |
% | dofIm = nodi interfaccia master |
% | dofIs = nodi interfaccia slave  |
% dof Lagrange multiplier not needed

% build an easy to access dof map
dofIm = DofMap.getCompDoF(nodesMaster);
dofM = DofMap.getCompDoF((1:masterMesh.nNodes)');
dofM = dofM(~ismember(dofM,dofIm));
dofIs = DofMap.getCompDoF(nodesSlave);
dofS = DofMap.getCompDoF((1:slaveMesh.nNodes)');
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
assert(nM==nS,'Meshes are not conforming!')
HD = zeros(2*nS,2*nS);
HM = zeros(2*nS,2*nM);
HD(2*(1:nS)'-1,2*(1:nS)'-1) = HDtmp;
HD(2*(1:nS)',2*(1:nS)') = HDtmp;
HM(2*(1:nS)'-1,2*(1:nM)'-1) = HMtmp;
HM(2*(1:nS)',2*(1:nM)') = HMtmp;

% complete saddle point matrix
K = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofIs));
   zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofIs)) ;
   KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -HM';
   zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, HD';
   zeros(length(dofIs), length(dofM)), zeros(length(dofIs), length(dofS)), -HM, HD,  zeros(length(dofIs), length(dofIs))
   ];

% creating global dof array and initializing forcing vector
listDofsSP = [dofM;dofS;dofIm; dofIs; dofIs];
nDofs = numel(listDofsSP);
f = zeros(nDofs,1);

% ------------------- APPLY BCS -------------------------------

if patch == 1
   % get Loaded dofs on top edge
   [nodesLoad,lInf] = getAreaInf(slaveMesh,2);
   % special treatment of extreme points (having force 1/2)
   loadDoF = dof.getDoF(nodesLoad,'slave',2,'y');
   f(loadDoF) = Fy*lInf;
end

%-------------------------LATERAL LOAD BCS ------------------------
if patch == 2 || patch == 3
   [nodesLoad,lInf] = getAreaInf(slaveMesh,3);
   loadDoF = dof.getDoF(nodesLoad,'slave',2,'x');
   f(loadDoF) = Fx*lInf;
end


%------------------- BOTTOM FIXED BCS -----------------------------
% get fixed dofs on bottom edge
% y bottom constraint
dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
if patch == 1 || patch == 2
   dirBotDoF = dof.getDoF(dirNod,'master',2);
elseif patch == 3
   dirBotDoF = dof.getDoF(dirNod,'master',2,'y');
end
[K,f] = applyDir(dirBotDoF, zeros(length(dirBotDoF),1), K, f);
% 
% 
% dirParas = dof.getDoF(123,'master',2,'y');
% [K,f] = applyDir(dirParas, 0, K, f);


%------------------- LATERAL FIXED BCS -----------------------------
if patch == 3
   dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 4,:));
   % remove interface slave node for unstabilized version
   if strcmp(stab,'unstable')
      [~,id] = min(slaveMesh.coordinates(dirNod,2));
      dirNod(id) = [];
   end
   dirDoF = dof.getDoF(dirNod,'slave',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
   %
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
   dirDoF = dof.getDoF(dirNod,'master',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
end



% remove lagrange multipliers belonging to the interface (then use constant
% interpolation and extend the adjacent value)
% get boundary nodes and adjacent node id
if strcmp(stab,'stable')
    [K,f] = handleEndPoints(slaveMesh,nodesSlave,dofM,dofS,dofIm,dofIs,K,f);
end

% solve linear system
u = K\f;

if strcmp(stab,'stable')
    u = addSolToEndPoints(slaveMesh,nodesSlave,dofM,dofS,dofIm,dofIs,u);
end

n = length(dofM)+length(dofS);
l = length(dofM)+length(dofS)+length(dofIm);
% master side
u_master = u(n+1:n+length(dofIm));
% get slave side solution
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
plotSolution(slaveMesh,'topConf',u_top,stressSlave);
plotSolution(masterMesh,'bottomConf',u_bottom,stressMaster);
%

% plot lagrange multipliers in matlab
% compute multipliers
[~,id] = sort(slaveMesh.coordinates(mortar.nodesSlave,1));
mult_x = mult(2*id-1);
mult_y = mult(2*id);
u_slave = [u_slave(2*id-1) u_slave(2*id)];
% fid = fopen('mult_fine.dat','w');
% for i=1:numel(mult_x)
%     fprintf(fid,'%1.9e %1.9e \n',mult_x(i),mult_y(i));
% end
[~,id] = sort(masterMesh.coordinates(mortar.nodesMaster,1));
u_master = [u_master(2*id-1) u_master(2*id)];
end

