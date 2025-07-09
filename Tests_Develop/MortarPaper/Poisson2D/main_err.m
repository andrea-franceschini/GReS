clc
clear
close all
% Mortar solution of a simple Poisson problem on a unit square domain
% exact solution is: u_ex = 16xy(1-x)(1-y);
% the rhs is f = 32x(1-x) + 32y(1-y);

% flagPlot: true --> results to Paraview in this directory
% fPlot = true;

% selecting solution method
% COND --> condansated approach
% SP --> saddle point matrix
sol_scheme = 'COND';

% IMPORT MESHES
masterMesh = Mesh();
slaveMesh = Mesh();

% Import meshes (no convergence profiles)
% topMesh.importGMSHmesh('Mesh/TopBlock_curve.msh');
% bottomMesh.importGMSHmesh('Mesh/BottomBlock_curve.msh');

% file Names for input meshes




% % generate meshes once forever
% for i = 1:nref
%   % write mesh to file
%   % set refinement
%   Nt = 1.5*(2*2^(i-1));
%   Nb = (2*2^(i-1));
%
%   fnameBot = strcat('bot_',num2str(i),'_',elem_type);
%   command = "python Mesh/scripts/domain_bot.py "  + fnameBot...
%     + " " + num2str(Nb) + " " + elem_type;
%   system(command);
%
%   fnameTop = strcat('top_',num2str(i),'_',elem_type);
%   command = "python Mesh/scripts/domain_top.py "  + fnameTop...
%     + " " + num2str(Nt) + " " + elem_type;
%   system(command);
%
% end


% selecting integration approach
integration_type = "SB";  % SB, RBF, EB
nInt = 6;
nGP = 4;

mult_type = 'dual';

elem_type = 'hexa';

nGrids = 3;
brokenL2 = zeros(nGrids,1);
brokenH1 = zeros(nGrids,1);
h = zeros(nGrids,1);

Nt = 13;
Nb = 5;


fnameBot = strcat('bot_err','_',elem_type);
command = "python Mesh/scripts/domain_bot.py "  + fnameBot...
  + " " + num2str(5) + " " + num2str(3) + " " + elem_type;
system(command);

fnameTop = strcat('top_err','_',elem_type);
command = "python Mesh/scripts/domain_top.py "  + fnameTop...
  + " " + num2str(13) + " " + num2str(7) + " " + elem_type;
system(command);

% Import the mesh data into the Mesh object
fnameBot = fullfile('Mesh','meshes',strcat('bot_err_',elem_type,'.vtk'));
fnameTop =  fullfile('Mesh','meshes',strcat('top_err_',elem_type,'.vtk'));

masterMesh.importMesh(fnameTop);
slaveMesh.importMesh(fnameBot);

% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh,nGP);
elemsSlave = Elements(slaveMesh,nGP);

% computing Stiffness matrix on top and bottom domainall
[KMaster, aNmaster] = stiffPoisson(masterMesh, elemsMaster);
[KSlave, aNslave] = stiffPoisson(slaveMesh, elemsSlave);

% get id of nodes belonging to master and slave interfaces
nodesMaster = unique(masterMesh.edges(masterMesh.edgeTag == 1,:));
nodesSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 1,:));

degree = 1;
if strcmp(elem_type,'hexa27')
  degree = 2;
end

% compute mortar operator
mortar = Mortar2D(degree,masterMesh,1,slaveMesh,1);

[D,M] = mortar.computeMortarSegmentBased(nGP,mult_type);

if strcmp(sol_scheme,'SP')
  % remove multipliers
  D(3,:) = D(3,:) + D(1,:);
  D(end,:) = D(end,:) + D(2,:);
  M(3,:) = M(3,:) + M(1,:);
  M(end,:) = M(end,:) + M(2,:);
  D(1,:) = []; D(1,:) = [];
  M(1,:) = []; M(1,:) = [];
end

% mortar operator
if strcmp(mult_type,'standard')
  E = D\M;
elseif strcmp(mult_type,'dual')
  E = (1./diag(D)).*M;
else
  error('wrong multiplier type string')
end


dofIm = nodesMaster;
dofM = (1:masterMesh.nNodes)';
dofM = dofM(~ismember(dofM,dofIm));
dofIs = nodesSlave;
dofS = (1:slaveMesh.nNodes)';
dofS = dofS(~ismember(dofS,dofIs));
Kmm = KMaster(dofM,dofM);
KmIm = KMaster(dofM,dofIm);
Kss = KSlave(dofS,dofS);
KsIs = KSlave(dofS, dofIs);
KImIm = KMaster(dofIm,dofIm);
KIsIs = KSlave(dofIs,dofIs);

%     hM = h(i);
%     H = hM*mortar.computePressureJumpMat();


% compute forcing vector (from analytical solution)

% get dofs coordinates
xM = [masterMesh.coordinates(dofM,1);
  masterMesh.coordinates(dofIm,1)];
yM = [masterMesh.coordinates(dofM,2);
  masterMesh.coordinates(dofIm,2)];
xS = [slaveMesh.coordinates(dofS,1);
  slaveMesh.coordinates(dofIs,1)];
yS = [slaveMesh.coordinates(dofS,2);
  slaveMesh.coordinates(dofIs,2)];
% compute forcing term (LM rows are set to 0)
x = [xM;xS];
y = [yM;yS];
fAn = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
%fAn = @(x,y) 32*(x.*(1-x)+y.*(1-y));
f = fAn(x,y);
% compute forcing vector and areanod
f = f.*[aNmaster(dofM);aNmaster(dofIm);aNslave(dofS);aNslave(dofIs)];
% reorder system rhs
f = [f(1:length(dofM));
  f(length(dofM)+length(dofIm)+1:length(dofM)+length(dofIm)+length(dofS));
  f(length(dofM)+1:length(dofM)+length(dofIm));
  f(length(dofM)+length(dofIm)+length(dofS)+1:end)];

if strcmp(sol_scheme, 'SP')
  % complete saddle point matrix
  K = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofIs)-2);
    zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofIs)-2) ;
    KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M';
    zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D';
    zeros(length(dofIs)-2, length(dofM)), zeros(length(dofIs)-2, length(dofS)), -M, D,  zeros(length(dofIs)-2,length(dofIs)-2)];

  f = [f; zeros(length(dofIs)-2,1)];
  listDofs = [dofM;dofS;dofIm; dofIs; dofIs];
else
  K = [Kmm, zeros(length(dofM),length(dofS)), KmIm;
    zeros(length(dofS),length(dofM)), Kss, KsIs*E;
    KmIm', E'*KsIs', KImIm+E'*KIsIs*E];
  listDofs = [dofM;dofS;dofIm];
  f(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm)) = ...
    f(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm)) + E'*f(length(dofM)+length(dofS)+length(dofIm)+1:end);
  f = f(1:length(dofM)+length(dofS)+length(dofIm));
end


% ------------------------------ APPLY BCS -------------------------------
% homogeneous Dirichlet BCs
% master domain
% get nodes from master domain (not in the interface)
nodesLatMaster = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
% extract nodes not belonging to the interface
nodesLatMaster = nodesLatMaster(~ismember(nodesLatMaster, nodesMaster));
% get corresponding DoFs in the linear system
dofLatMaster = nodesLatMaster;
dofLatMaster = getGlobalDofs(dofLatMaster,dofM,dofS,dofIm,dofIs,'master');
% Apply penalty method
[K,f] = applyDir(dofLatMaster, zeros(length(dofLatMaster),1), K, f);

% slave domain
nodesLatSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
% extract nodes not belonging to the interface
nodesLatSlave = nodesLatSlave(~ismember(nodesLatSlave, nodesSlave));
% get corresponding DoFs in the linear system
dofLatSlave = nodesLatSlave; % x direction is fixed
dofLatSlave = getGlobalDofs(dofLatSlave,dofM,dofS,dofIm,dofIs,'slave');
% Apply Dirichlet BCs
[K,f] = applyDir(dofLatSlave, zeros(length(dofLatSlave),1), K, f);

% master interface
nodesLatMaster = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
nodesLatInt = nodesLatMaster(ismember(nodesLatMaster, nodesMaster));
dofLatInt = nodesLatInt;
dofLatInt = getGlobalDofs(dofLatInt,dofM,dofS,dofIm,dofIs,'interfaceMaster');
% Apply Dirichlet BCs
[K,f] = applyDir(dofLatInt, zeros(length(dofLatInt),1), K, f);
%
% -------------------------------- SOLVE SYSTEM --------------------------

% solve linear system
u = K\f;

u_master = zeros(masterMesh.nNodes,1);
u_slave = zeros(slaveMesh.nNodes,1);
%collect displacement of master domain and slave domain, according to user assignment;
u_master(dofM) = u(1:length(dofM));
u_master(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
l = length(dofM)+length(dofS)+length(dofIm);
if strcmp(sol_scheme,'SP')
  u_s = u(l+1:l+length(dofIs));
else
  u_s = E*u(length(dofM)+length(dofS)+1:l);
end


u_slave(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
u_slave(dofIs) = u_s;


%% Plotting solutions
u_anal = @(x,y) sin(pi*x).*sin(pi*y);
dux =  @(x,y) pi*cos(pi*x).*sin(pi*y);
duy = @(x,y) pi*sin(pi*x).*cos(pi*y);

% u_anal = @(x,y) 16*x.*y.*(1-x).*(1-y);
% dux =  @(x,y) 16*y - 16*y.^2 - 32*x.*y + 32*x.*y.^2;
% duy = @(x,y) 16*x - 16*x.^2 - 32*x.*y + 32*y.*x.^2;

c = masterMesh.coordinates;
u_anal_master = u_anal(c(:,1),c(:,2));
plotFunction(masterMesh,'OUT_errMaster',abs(u_master - u_anal_master));
c = slaveMesh.coordinates;
u_anal_slave = u_anal(c(:,1),c(:,2));
plotFunction(slaveMesh,'OUT_errSlave',abs(u_slave - u_anal_slave));




