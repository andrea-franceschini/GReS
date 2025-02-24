function [mult_x,mult_y,u_slave,u_master,x] = RunNonConfConstant(Fx,Fy,DmatS,DmatM,patch,nXs,nYs,nXm,nYm,Em,Es)
% Run a 2 non conforming block patch test with constant multipliers

% plots solution to paraview
% output: multipliers on the interface

% gauss class for 1D element integration
nG = 16;

% Set the input file name
gaussQuad = Gauss(12,2,2);


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
cellsSlave = find(slaveMesh.edgeTag == 1);
cellsMaster = find(masterMesh.edgeTag == 1);


dof = DofMap(masterMesh,nodesMaster,slaveMesh,nodesSlave);
% % remove extreme nodes in contact with the boundary
mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
% compute mortar operator
[D,M] = mortar.computeMortarConstant(nG,6);
D = expandMat(D,2);
M = expandMat(M,2);
hS = 1/numel(cellsSlave); % mesh size
hM = 1/numel(cellsMaster); % mesh size
E = 0.5*(Em+Es);
H = mortar.computePressureJumpMat();
H = expandMat(H,2);


% build an easy to access dof map
dofIm = DofMap.getCompDoF(nodesMaster);
dofM = DofMap.getCompDoF((1:masterMesh.nNodes)');
dofM = dofM(~ismember(dofM,dofIm));
dofIs = DofMap.getCompDoF(nodesSlave);
dofMult = DofMap.getCompDoF(cellsSlave);
dofS = DofMap.getCompDoF((1:slaveMesh.nNodes)');
dofS = dofS(~ismember(dofS,dofIs));
Kmm = KMaster(dofM,dofM);
KmIm = KMaster(dofM,dofIm);
Kss = KSlave(dofS,dofS);
KsIs = KSlave(dofS, dofIs);
KImIm = KMaster(dofIm,dofIm);
KIsIs = KSlave(dofIs,dofIs);

stab = 0.5*hS/E;
H = stab*H;
%H = zeros(numel(dofMult),numel(dofMult));

% pressure-jump stabilization matrix ( regular matrix with homogeneous
% material)
r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofMult))];
r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofMult))] ;
r3 = [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M'];
r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D'];
r5 = [zeros(length(dofMult), length(dofM)), zeros(length(dofMult), length(dofS)), -M, D, -H];
K = [r1;r2;r3;r4;r5];



% complete saddle point matrix
% creating global dof array and initializing forcing vector
listDofsSP = [dofM;dofS;dofIm;dofIs;dofMult];
f = zeros(length(listDofsSP),1);


% ------------------- APPLY BCS -------------------------------

%------------------- TOP LOAD BCS -----------------------------

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

%------------------- LATERAL FIXED BCS -----------------------------
if patch == 3
   dirNod = unique(slaveMesh.edges(slaveMesh.edgeTag == 4,:));
   % remove interface slave node for unstabilized version
   dirDoF = dof.getDoF(dirNod,'slave',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
   %
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
   dirDoF = dof.getDoF(dirNod,'master',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
end

% masterDof = dof.getDoF(nodesMaster(3:end),'master');
% slaveDof = dof.getDoF(nodesSlave(3:end),'slave');
% lagDof = dof.getDoF(nodesSlave(3:end),'lagrange');


% solve linear system
u = K\f;
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
plotSolution(slaveMesh,'topConstant',u_top,stressSlave);
plotSolution(masterMesh,'botConstant',u_bottom,stressMaster);
%
% get stress on slave interface
% cellInterf = any(ismember(slaveMesh.surfaces,nodesSlave),2);
% s_x = stressSlave(cellInterf,3);
% s_y = stressSlave(cellInterf,2);


% plot lagrange multipliers in matlab
% compute multipliers
[~,id] = sort(masterMesh.coordinates(mortar.nodesMaster,1));
u_master = [u_master(2*id-1) u_master(2*id)];


[~,id] = sort(slaveMesh.coordinates(mortar.nodesSlave,1));
u_slave = [u_slave(2*id-1) u_slave(2*id)];

% reorder edges on the slave side
c = [slaveMesh.coordinates(slaveMesh.edges(cellsSlave,1),1), slaveMesh.coordinates(slaveMesh.edges(cellsSlave,2),1)]; 
x = 0.5*sum(c,2);
[x,id] = sort(x);
mult_x = mult(2*id-1);
mult_y = mult(2*id);

end

