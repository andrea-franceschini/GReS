function [mult_x,mult_y,ux_slave,uy_slave,coordInt] = RunConfPatch2(Fx,Fy,Dmat,patch,stab,nX,nY,nTip)
% Run a 2 non conforming block patch test with standard or dual lagrange
% multipliers

% plots solution to paraview
% output: multipliers on the interface

% DEFINE MODEL

% IMPORT MESHES
if nTip==0
    masterMesh = getMesh('Mesh_flat/bottomBlock.geo','bottom',nX,nY);
    slaveMesh = getMesh('Mesh_flat/topBlock.geo','top',nX,nY);
else
    masterMesh = getMesh('Mesh_flat/bottomBlockTip.geo','bottom',nX,nY,nX,nTip);
    slaveMesh = getMesh('Mesh_flat/topBlockTip.geo','top',nX,nY,nX,nTip);
end

% Set the input file name
gaussQuad = Gauss(12,2,2);


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
listDofsSP = [dofM';dofS';dofIm'; dofIs'; dofIs'];
nDofs = numel(listDofsSP);
f = zeros(nDofs,1);

% ------------------- APPLY BCS -------------------------------

%------------------- TOP LOAD BCS -----------------------------
% get Loaded dofs on top edge
switch patch
    case 1
        nodesLoad = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
        % special treatment of infextreme points (having force 1/2)
        n_ext = nodesLoad(ismember(slaveMesh.coordinates(nodesLoad,1),[0; 1]));
        loadDoFext = 2*n_ext; % loading y component
        loadDoFY = 2*nodesLoad(~ismember(nodesLoad, n_ext)); % loading y component
        loadDoFext = getGlobalDofs(loadDoFext, dofM,dofS,dofIm,dofIs, 'slave');
        loadDoFY = getGlobalDofs(loadDoFY, dofM,dofS,dofIm,dofIs, 'slave');
        f(loadDoFext) = Fy*(0.5/(length(nodesLoad)-1));
        f(loadDoFY) = Fy*(1/(length(nodesLoad)-1));
end

%-------------------------LATERAL LOAD BCS ------------------------
switch patch
    case {2,3}
        nodesLoad = find(abs(slaveMesh.coordinates(:,1)-0)<1e-3);
        % get top node (it has half of the entities influence)
        [yC,id] = sort(slaveMesh.coordinates(nodesLoad,2),'ascend');
        l1 = abs(yC(end)-yC(end-1));
        l2 = abs(yC(2)-yC(1));
        n1 = nodesLoad(id(end));
        f(getGlobalDofs(2*n1-1,dofM,dofS,dofIm,dofIs,'slave')) = Fx*l1/2;
        nIn = nodesLoad(id(3:end-1));
        f(getGlobalDofs(2*nIn-1,dofM,dofS,dofIm,dofIs,'slave')) = Fx*l1;
        n2 = nodesLoad(id(2));
        f(getGlobalDofs(2*n2-1,dofM,dofS,dofIm,dofIs,'slave')) = Fx*(l1/2+l2/2);
        n3 = nodesLoad(id(1));
        f(getGlobalDofs(2*n3-1,dofM,dofS,dofIm,dofIs,'interfaceSlave')) = Fx*l2/2;
        %
        % % LATERAL LOAD ALSO TO MASTER SIDE (BOTTOM)
        %get top node (it has half of the entities influence)
%         nodesLoad = find(abs(masterMesh.coordinates(:,1)-0)<1e-3);
%         [yC,id] = sort(masterMesh.coordinates(nodesLoad,2),'descend');
%         l = abs(yC(end)-yC(end-1));
%         n1 = nodesLoad(id(1));
%         f(getGlobalDofs(2*n1-1,dofM,dofS,dofIm,dofIs,'interfaceMaster')) = Fx*l/2;
%         nIn = nodesLoad(id(2:end-1));
%         f(getGlobalDofs(2*nIn-1,dofM,dofS,dofIm,dofIs,'master')) = Fx*l;
%         n2 = nodesLoad(id(end));
%         f(getGlobalDofs(2*n2-1,dofM,dofS,dofIm,dofIs,'master')) = Fx*l/2;
end

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
switch patch
    case 3
        % get nodes on right edge of slave domain
        nodesFixX = find(abs(slaveMesh.coordinates(:,1)-1)<1e-3);
        dirXDoF = getGlobalDofs(2*nodesFixX-1,dofM,dofS,dofIm,dofIs,'slave');
        id = abs(slaveMesh.coordinates(nodesFixX,2)-1)<1e-3;
        dirXIntDoF = getGlobalDofs(2*nodesFixX(id)-1,dofM,dofS,dofIm,dofIs,'interfaceSlave');
        [K,f] = applyDir(dirXDoF, zeros(length(dirXDoF),1), K, f);
        if strcmp(stab,'stable')
            [K,f] = applyDir(dirXIntDoF, zeros(length(dirXIntDoF),1), K, f);
        end

        % get nodes on right edge of master domain
        nodesFixX = find(abs(masterMesh.coordinates(:,1)-1)<1e-3);
        dirXDoF = getGlobalDofs(2*nodesFixX-1,dofM,dofS,dofIm,dofIs,'master');
        id = abs(masterMesh.coordinates(nodesFixX,2)-1)<1e-3;
        dirXIntDoF = getGlobalDofs(2*nodesFixX(id)-1,dofM,dofS,dofIm,dofIs,'interfaceMaster');
        [K,f] = applyDir(dirXDoF, zeros(length(dirXDoF),1), K, f);
        [K,f] = applyDir(dirXIntDoF, zeros(length(dirXIntDoF),1), K, f);
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
stressMaster = computeStress(masterMesh,elemsMaster,Dmat,u_bottom,gaussQuad);
stressSlave = computeStress(slaveMesh,elemsSlave,Dmat,u_top,gaussQuad);
%
plotSolution(slaveMesh,'topConf',u_top,stressSlave);
plotSolution(masterMesh,'bottomConf',u_bottom,stressMaster);
%


% plot lagrange multipliers in matlab
% compute multipliers
[coordInt,id] = sort(slaveMesh.coordinates(mortar.nodesSlave,1));
mult_x = mult(2*id-1);
mult_y = mult(2*id);
ux_slave = u_slave(2*id-1);
uy_slave = u_slave(2*id);
end

