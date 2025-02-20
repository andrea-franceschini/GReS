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
listDofsSP = [dofM';dofS';dofIm'; dofIs'; dofIm'; dofIs'];
f = zeros(length(listDofsSP),1);


% ------------------- APPLY BCS -------------------------------

%------------------- TOP LOAD BCS -----------------------------

% get Loaded dofs on top edge
nodesLoad = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
% special treatment of infextreme points (having force 1/2)
n_ext = nodesLoad(ismember(slaveMesh.coordinates(nodesLoad,1),[0; 1]));
loadDoFext = 2*n_ext; % loading y component
loadDoFY = 2*nodesLoad(~ismember(nodesLoad, n_ext)); % loading y component
loadDoFext = getGlobalDofs(loadDoFext, dofM,dofS,dofIm,dofIs, 'slave');
loadDoFY = getGlobalDofs(loadDoFY, dofM,dofS,dofIm,dofIs, 'slave');
f(loadDoFext) = Fy*(0.5/(length(nodesLoad)-1));
f(loadDoFY) = Fy*(1/(length(nodesLoad)-1));

%-------------------------LATERAL LOAD BCS ------------------------
nodesLoad = find(abs(slaveMesh.coordinates(:,1)-0)<1e-3);
% get top node (it has half of the entities influence)
[yC,id] = sort(slaveMesh.coordinates(nodesLoad,2),'ascend');
l1 = abs(yC(end)-yC(end-1));
l2 = abs(yC(2)-yC(1));
n1 = nodesLoad(id(end));
f(getGlobalDofs(2*n1-1,dofM,dofS,dofIm,dofIs,'slave')) = Fx*l1/2;
nIn = nodesLoad(id(3:end-1));
f(getGlobalDofs(2*nIn-1,dofM,dofS,dofIm,dofIs,'slave')) = Fx*l1;
c2 = nodesLoad(id(2));
f(getGlobalDofs(2*c2-1,dofM,dofS,dofIm,dofIs,'slave')) = Fx*(l1/2+l2/2);
n3 = nodesLoad(id(1));
f(getGlobalDofs(2*n3-1,dofM,dofS,dofIm,dofIs,'interfaceSlave')) = Fx*l2/2;
%
% % LATERAL LOAD ALSO TO MASTER SIDE (BOTTOM)
%get top node (it has half of the entities influence)

nodesLoad = find(abs(masterMesh.coordinates(:,1)-0)<1e-3);
[yC,id] = sort(masterMesh.coordinates(nodesLoad,2),'descend');
l = abs(yC(end)-yC(end-1));
n1 = nodesLoad(id(1));
%f(getGlobalDofs(2*n1-1,dofM,dofS,dofIm,dofIs,'interfaceMaster')) = Fx*l/2;
nIn = nodesLoad(id(2:end-1));
%f(getGlobalDofs(2*nIn-1,dofM,dofS,dofIm,dofIs,'master')) = Fx*l;
c2 = nodesLoad(id(end));
%f(getGlobalDofs(2*n2-1,dofM,dofS,dofIm,dofIs,'master')) = Fx*l/2;

%------------------- BOTTOM FIXED BCS -----------------------------
% get fixed dofs on bottom edge
% y bottom constraint
dirNod = find(abs(masterMesh.coordinates(:,2)-0)<1e-3);
dirBotDoFY = getGlobalDofs(2*dirNod,dofM,dofS,dofIm,dofIs,'master');
dirBotDoFX = getGlobalDofs(2*dirNod-1,dofM,dofS,dofIm,dofIs,'master');
[K,f] = applyDir(dirBotDoFY, zeros(length(dirBotDoFY),1), K, f);
switch patch
    case {1,2}
        [K,f] = applyDir(dirBotDoFX, zeros(length(dirBotDoFX),1), K, f);
end


% -------------------LATERAL CONSTRAINT BCS-----------------------
% get nodes on right edge of slave domain
switch patch
    case 3
        nodesFixX = find(abs(slaveMesh.coordinates(:,1)-1)<1e-3);
        dirXDoF = getGlobalDofs(2*nodesFixX-1,dofM,dofS,dofIm,dofIs,'slave');
        id = abs(slaveMesh.coordinates(nodesFixX,2)-1)<1e-3;
        dirXIntDoF = getGlobalDofs(2*nodesFixX(id)-1,dofM,dofS,dofIm,dofIs,'interfaceSlave');
        [K,f] = applyDir(dirXDoF, zeros(length(dirXDoF),1), K, f);
        % get nodes on right edge of master domain
        nodesFixX = find(abs(masterMesh.coordinates(:,1)-1)<1e-3);
        dirXDoF = getGlobalDofs(2*nodesFixX-1,dofM,dofS,dofIm,dofIs,'master');
        id = abs(masterMesh.coordinates(nodesFixX,2)-1)<1e-3;
        dirXIntDoF = getGlobalDofs(2*nodesFixX(id)-1,dofM,dofS,dofIm,dofIs,'interfaceMaster');
        [K,f] = applyDir(dirXDoF, zeros(length(dirXDoF),1), K, f);
        [K,f] = applyDir(dirXIntDoF, zeros(length(dirXIntDoF),1), K, f);
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

