% 2D linear mechanics patch test

% INPUT
nG = 15;
nInt = 4;
gaussQuad = Gauss(12,2,2);

patch = 3;       % chose patch test case
switch patch
    case 1
        Fx = 0; %[kPa]
        Fy = -10;
        nu = 0;
    case {2,3}
        Fx = 10; %[kPa]
        Fy = 0;
        nu = 0.25;
end

% Elastic properties and constitutive tensor
% Mortar side
E1 = 100000;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E1/((1+nu)*(1-2*nu)))*D;
DmatM = D;
% Slave side
E2 = 10000;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E2/((1+nu)*(1-2*nu)))*D;
DmatS = D;


nel = 60;   % number of elements on the slave interface
rat = 0.3;  % numb master / numb slave elems

%fig = figure('Visible', 'off');
fig.Position = [100,100,800,600];

t = tiledlayout(1,2);
nexttile(1)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t_x$', 'Interpreter', 'latex')
if patch == 1
   ylim([-0.01 0.01])
end
hold on
nexttile(2)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$t_y$', 'Interpreter', 'latex')
if patch == 1
   ylim([-10.01 -9.99])
end
hold on


for scheme = ["DUAL","STANDARD","CONFORMING","P0"] % accomodate UNBIASED in the future!!!
nXs = nel+1;
nYs = round(1*nXs);
nXm = round(nel*rat+1);
if strcmp(scheme,'CONFORMING')
   nXm = nXs;
end
nYm = round(1*nXm);
% Import the mesh data into the Mesh object
masterMesh = getMesh('Mesh/bottomBlock.geo','bottom',nXm,nYm);
slaveMesh = getMesh('Mesh/topBlock.geo','top',nXs,nYs);
% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh,gaussQuad);
elemsSlave = Elements(slaveMesh,gaussQuad);

Emaster = E1*ones(masterMesh.nSurfaces,1); 
Eslave = E2*ones(slaveMesh.nSurfaces,1); 

% Mesh size
h = 1/nXs;

mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);

switch scheme
   case 'DUAL'
      [D,M] = mortar.computeMortarSegmentBased(nG,'dual');
      leg = 'Dual';
   case 'STANDARD'
      [D,M] = mortar.computeMortarSegmentBased(nG,'standard');
      leg = 'Standard';
   case 'P0'
      [D,M] = mortar.computeMortarConstant(nG,nInt);
      leg = 'P0';
   case 'CONFORMING'
      [D,M] = mortar.computeConfCrossGridMat();
      leg = 'Conforming';
end

nMult = 2*numel(mortar.nodesSlave);

% handle lagrange multipliers at the interface boundary
switch scheme
   case {'DUAL','STANDARD','CONFORMING'}
      E = D\M;
      D = handleEndPoints(D,mortar);
      M = handleEndPoints(M,mortar);
      E = expandMat(E,2);
      nMult = nMult - 4;
   case {'P0'}
      nMult = 2*mortar.nElSlave;
end

alpha = 10; % scalar empirical stabilization parameter for unbiased formulation
gamma = (alpha*h)/E1;

dof = DofMap(masterMesh,mortar.nodesMaster,slaveMesh,mortar.nodesSlave);
dofIm = DofMap.getCompDoF(mortar.nodesMaster,2);
dofM = DofMap.getCompDoF((1:masterMesh.nNodes)',2);
dofM = dofM(~ismember(dofM,dofIm));
dofIs = DofMap.getCompDoF(mortar.nodesSlave,2);
dofS = DofMap.getCompDoF((1:slaveMesh.nNodes)',2);
dofS = dofS(~ismember(dofS,dofIs));

KMaster = stiff(masterMesh, elemsMaster, DmatM, gaussQuad);
KSlave = stiff(slaveMesh, elemsSlave, DmatS, gaussQuad);
% get id of nodes belonging to master and slave interfaces
Kmm = KMaster(dofM,dofM);
KmIm = KMaster(dofM,dofIm);
Kss = KSlave(dofS,dofS);
KsIs = KSlave(dofS, dofIs);
KImIm = KMaster(dofIm,dofIm);
KIsIs = KSlave(dofIs,dofIs);

if strcmp(scheme,'P0')
   H1 = (h/E1)*mortar.computePressureJumpMat();
   H1 = full(expandMat(H1,2));
   H2 = mortar.computeStabilizationMatrix(Emaster,Eslave);
   H = H2;
   % h = vecnorm(H2,2,2);
   %H = diag(h)*H1;
else
   H = zeros(nMult,nMult);
end
D = expandMat(D,2);
M = expandMat(M,2);

%H = zeros(nMult,nMult);
f = zeros(size(KMaster,1)+size(KSlave,1),1);

switch scheme
   case 'UNBIASED'
      s = (hM);
      r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofIm)), zeros(length(dofM),length(dofIs))];
      r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofIm)), zeros(length(dofS),length(dofIs))] ;
      r3 =  [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), 0.5*D2', -0.5*M1'];
      r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, -0.5*M2', 0.5*D1'];
      r5 = [zeros(length(dofIm), length(dofM)), zeros(length(dofIm), length(dofS)), 0.5*D2, -0.5*M2,  s*D2, s*M2];
      r6 = [zeros(length(dofIs), length(dofM)), zeros(length(dofIs), length(dofS)), -0.5*M1, 0.5*D1,  s*M1, s*D1];
      K = [r1;r2;r3;r4;r5;r6];
      f = [f; zeros(length(dofIm)+length(dofIs),1)];
   otherwise
      r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),nMult)];
      r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),nMult)] ;
      r3 = [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M'];
      r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D'];
      r5 = [zeros(nMult, length(dofM)), zeros(nMult, length(dofS)), -M, D, -H];
      K = [r1;r2;r3;r4;r5];
      f = [f; zeros(nMult,1)];
end

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
   % [~,id] = min(slaveMesh.coordinates(dirNod,2));
   % dirNod(id) = [];
   dirDoF = dof.getDoF(dirNod,'slave',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
   %
   dirNod = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
   dirDoF = dof.getDoF(dirNod,'master',2,'x');
   [K,f] = applyDir(dirDoF, zeros(length(dirDoF),1), K, f);
end

% solve linear system
u = K\f;

% get multipliers on the interface
mult = u(end-nMult+1:end);



switch scheme
   case {'DUAL','CONFORMING','UNBIASED','STANDARD'}
      [x,id] = sort(slaveMesh.coordinates(mortar.nodesSlave(3:end),1));
   case 'P0'
      c = [slaveMesh.coordinates(mortar.slaveTopol(:,1),1), slaveMesh.coordinates(mortar.slaveTopol(:,2),1)];
      x = 0.5*sum(c,2);
      [x,id] = sort(x);
      %collect displacement of master domain and slave domain, according to user assignment;
      u_master = u(DofMap.getCompDoF(dof.nodeMapMaster));
      u_slave = u(DofMap.getCompDoF(dof.nodeMapSlave));
      plotSolution(slaveMesh,'slaveOut',u_slave);
      plotSolution(masterMesh,'masterOut',u_master);
      %
end

tx = mult(1:2:end);
ty = mult(2:2:end);

% PLOT MULTIPLIER
nexttile(1)
plot(x,tx(id),'s-','LineWidth',1,'DisplayName',leg)
hold on

nexttile(2)
plot(x,ty(id),'s-','LineWidth',1,'DisplayName',leg)
hold on

% naive projection of nodal multipliers
% switch scheme
%    case 'DUAL'
%       mortarSwap = Mortar2D(1,slaveMesh,1,masterMesh,1);
%       [Dsw, Msw] = mortarSwap.computeMortarSegmentBased(2,'dual');
%    case 'STANDARD'
%       mortarSwap = Mortar2D(1,slaveMesh,1,masterMesh,1);
%       [Dsw, Msw] = mortarSwap.computeMortarSegmentBased(2,'standard');
% end
% 
% switch scheme
%    case {'DUAL','STANDARD'}
%       Dsw = handleEndPoints(Dsw,mortarSwap);
%       Msw = handleEndPoints(Msw,mortarSwap);
%       Esw = Dsw\Msw;
%       Esw = expandMat(Esw,2);
%       mult_sw = [mult(1);mult(2);mult(end-1);mult(end-2);mult];
%       mult_proj = E*(Esw*mult_sw);
%       tx = mult_proj(1:2:end);
%       ty = mult_proj(2:2:end);
%       % plot projected multiplier
%       [x,id] = sort(masterMesh.coordinates(mortarSwap.nodesSlave,1));
%       nexttile(1)
%       plot(x,tx(id),'s-','LineWidth',1,'DisplayName',strcat('proj-',scheme))
%       nexttile(2)
%       plot(x,ty(id),'s-','LineWidth',1,'DisplayName',strcat('proj-',scheme))
% end

end


% finilize plot
nexttile(1)
legend show
set(legend, 'Interpreter', 'latex'); % Set legend to LaTeX
legend('Location', 'best');

nexttile(2)
legend show
set(legend, 'Interpreter', 'latex'); % Set legend to LaTeX
legend('Location', 'best');

fig.Visible = "on";
%exportgraphics(fig, strcat('Plots/Patch_',num2str(patch),'.pdf'), 'Resolution', 300); % High resolution


function mat = handleEndPoints(mat,mortar)
% return mortar matrix with removed multipliers at the tip
[~,i] = sort(mortar.slaveCoord(mortar.nodesSlave,1));
mat(i(2),:) = mat(i(2),:)+mat(i(1),:);
mat(i(end-1),:) = mat(i(end-1),:)+mat(i(end),:);
mat([i(1) i(end)],:) = [];
end