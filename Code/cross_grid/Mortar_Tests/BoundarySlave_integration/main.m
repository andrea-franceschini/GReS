clear; close all;

% Patch test for  boundary slave integration 

top = 'master';
nG = 20;
p = -0.5; % vertical pressure
E = 100; 
nu = 0;
integration = "P0"; % RBR,SB,P0

Dmat = zeros(3);
Dmat([1 5]) = 1-nu;
Dmat([2 4]) = nu;
Dmat(9) = 0.5*(1-2*nu);
Dmat = (E/((1+nu)*(1-2*nu)))*Dmat;

% Set the input file name
gaussQuad = Gauss(12,2,2);


% read mesh
[bottomMesh,topMesh] = deal(Mesh(),Mesh()); 
bottomMesh.importGMSHmesh('Mesh/botBlock.msh')
topMesh.importGMSHmesh('Mesh/topBlock.msh');

% Element class for further stiffness matrix computation
elemsTop = Elements(topMesh,gaussQuad);
elemsBot = Elements(bottomMesh,gaussQuad);


% assignign blocks to master and slave side
if strcmp(top,'master')
   masterMesh = topMesh;
   slaveMesh = bottomMesh;
   elemsMaster = elemsTop;
   elemsSlave = elemsBot;
   bottom = 'slave';
elseif strcmp(top,'slave')
   masterMesh = bottomMesh;
   slaveMesh = topMesh;
   elemsMaster = elemsBot;
   elemsSlave = elemsTop;
   bottom = 'master';
end

KMaster = stiff(masterMesh, elemsMaster, Dmat, gaussQuad);
KSlave = stiff(slaveMesh, elemsSlave, Dmat, gaussQuad);

mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
dof = DofMap(masterMesh,mortar.nodesMaster,slaveMesh,mortar.nodesSlave);


% compute mortar operator
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

% creating global dof array and initializing forcing vector
listDofsSP = [dofM;dofS;dofIm; dofIs;dofMult];
f = zeros(length(listDofsSP),1);


%------------------- TOP LOAD BCS -----------------------------

% get Loaded dofs on top edge
[nodesLoad,lInf] = getAreaInf(topMesh,2);
% special treatment of extreme points (having force 1/2)
loadDoF = dof.getDoF(nodesLoad,top,'y');
f(loadDoF) = p*lInf;

% load free nodes of bottom free face
% force manually computed
[nodesLoad,lInf] = getAreaInf(bottomMesh,1);
[~,id] = sort(bottomMesh.coordinates(nodesLoad,1));
nodesLoad = nodesLoad(id);
loadDoF = dof.getDoF(nodesLoad,bottom,'y');
% manually applying correct force to bottom nodes
f(loadDoF) = p*[0.1;0.12775;2.25e-3;2.25e-3;0.12775;0.1];



%------------------- BOTTOM FIXED BCS -----------------------------
% get fixed dofs on bottom edge
% y bottom constraint
dirNod = unique(bottomMesh.edges(bottomMesh.edgeTag == 2,:));
dirBotDoF = dof.getDoF(dirNod,bottom);
[K,f] = applyDir(dirBotDoF, zeros(length(dirBotDoF),1), K, f);

% solve linear system
u = K\f;

n = length(dofM)+length(dofS);
l = length(dofM)+length(dofS)+length(dofIm);
% get slave side solution
u_master = u(n+1:n+length(dofIm));
u_slave = u(l+1:l+length(dofIs));
mult = u(l+numel(dofIs)+1:end);

%collect displacement of master domain and slave domain, according to user assignment;
switch top
   case 'slave'
      u_top = zeros(slaveMesh.nNodes,1);
      u_bottom = zeros(masterMesh.nNodes,1);
      u_bottom(dofM) = u(1:length(dofM));
      u_bottom(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
      u_top(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
      u_top(dofIs) = u_slave;
   case 'master'
      u_bottom = zeros(slaveMesh.nNodes,1);
      u_top = zeros(masterMesh.nNodes,1);
      u_top(dofM) = u(1:length(dofM));
      u_top(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
      u_bottom(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
      u_bottom(dofIs) = u_slave;

end

% compute stresses
stressTop = computeStress(topMesh,elemsTop,Dmat,u_top,gaussQuad);
stressBottom = computeStress(bottomMesh,elemsBot,Dmat,u_bottom,gaussQuad);
%
plotSolution(bottomMesh,'bottom',u_bottom,stressBottom);
plotSolution(topMesh,'top',u_top,stressTop);
%
switch integration
   case {'RBF','SB'}
      [x,id] = sort(slaveMesh.coordinates(mortar.nodesSlave,1));
   case 'P0'
      c = [slaveMesh.coordinates(mortar.slaveTopol(:,1),1), slaveMesh.coordinates(mortar.slaveTopol(:,1),2)];
      x = 0.5*sum(c,2);
      [x,id] = sort(x);
end
mult_x = mult(2*id-1);
mult_y = mult(2*id);

m = -0.5;
switch top
   case 'master'
      m = 0.5;
end

figure(1)
plot([min(x) max(x)],[m,m],'k-','LineWidth',1.2)
hold on
plot(x,mult_y,'b-s')





