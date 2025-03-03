function [mult_x,mult_y,ux_slave,uy_slave] = RunPuso2020(Fx,Fy,DmatS,DmatM,gamma,patch,nXs,nYs,nXm,nYm)
% Run a 2 non conforming block patch test according to Puso 2020
% unbiased stabilized formulation.
% stick to standard multiplier case
% gamma: stabilization parameter

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

dof = DofMap(masterMesh,nodesMaster,slaveMesh,nodesSlave);
mortar1 = Mortar2D(1,masterMesh,1,slaveMesh,1);
mortar2 = Mortar2D(1,slaveMesh,1,masterMesh,1);
% compute mortar operator

[D1,M1] = mortar1.computeMortarRBF(nG,4,'gauss','standard');
[D2,M2] = mortar2.computeMortarRBF(nG,4,'gauss','standard');

h = 1/(0.5*(numel(nodesSlave)+numel(nodesMaster)));

% reordering the matrix of the system
%
% | dofM = nodi interni master      |
% | dofS = nodi interni slave       |
% | dofIm = nodi interfaccia master |
% | dofIs = nodi interfaccia slave  |
% | moltiplicatori master           |
% | moltiplicatori slave            |

% build an easy to access dof map
dofIm = DofMap.getCompDoF(nodesMaster');
dofM = DofMap.getCompDoF(1:masterMesh.nNodes);
dofM = dofM(~ismember(dofM,dofIm));
dofIs = DofMap.getCompDoF(nodesSlave');
dofS = DofMap.getCompDoF(1:slaveMesh.nNodes);
dofS = dofS(~ismember(dofS,dofIs));
Kmm = KMaster(dofM,dofM);
KmIm = KMaster(dofM,dofIm);
Kss = KSlave(dofS,dofS);
KsIs = KSlave(dofS, dofIs);
KImIm = KMaster(dofIm,dofIm);
KIsIs = KSlave(dofIs,dofIs);
% expanding the mortar operator also to the y direction
D1 = expandMat(D1,2);
D2 = expandMat(D2,2);
M1 = expandMat(M1,2);
M2 = expandMat(M2,2);

s = 0.5*gamma*h;
% complete saddle point matrix
r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofIm)), zeros(length(dofM),length(dofIs))];
r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofIm)), zeros(length(dofS),length(dofIs))] ;
r3 =  [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), 0.5*D2', -0.5*M1'];
r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, -0.5*M2', 0.5*D1'];
r5 = [zeros(length(dofIm), length(dofM)), zeros(length(dofIm), length(dofS)), 0.5*D2, -0.5*M2,  s*D2, s*M2];
r6 = [zeros(length(dofIs), length(dofM)), zeros(length(dofIs), length(dofS)), -0.5*M1, 0.5*D1,  s*M1, s*D1];
K = [r1;r2;r3;r4;r5;r6];
% creating global dof array and initializing forcing vector
listDofsSP = [dofM;dofS;dofIm;dofIs;dofIm;dofIs];
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
   [~,id] = min(slaveMesh.coordinates(dirNod,2));
   dirNod(id) = [];
   dirDoF = dof.getDoF(dirNod,'slave',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
   %
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
   dirDoF = dof.getDoF(dirNod,'master',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
end


u = K\f;


n = length(dofM)+length(dofS);
l = length(dofM)+length(dofS)+length(dofIm);
% get slave side solution

u_slave = u(l+1:l+length(dofIs));
mult = u(end-numel(u_slave)+1:end);


%mortarSwap.modifyEndPoints();

% plot lagrange multipliers in matlab
% compute multipliers
[~,id] = sort(slaveMesh.coordinates(nodesSlave,1));
mult_x = mult(2*id-1);
mult_y = mult(2*id);
ux_slave = u_slave(2*id-1);
uy_slave = u_slave(2*id);


end

