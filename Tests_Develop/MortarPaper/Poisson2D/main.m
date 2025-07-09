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

mult_type = 'dual'; 

% IMPORT MESHES
masterMesh = Mesh();
slaveMesh = Mesh();

% Import meshes (no convergence profiles)
% topMesh.importGMSHmesh('Mesh/TopBlock_curve.msh');
% bottomMesh.importGMSHmesh('Mesh/BottomBlock_curve.msh');

% file Names for input meshes
elem_type = 'hexa27';
switch elem_type
  case {'hexa','tetra'}
    degree = 1;
  case 'hexa27'
    degree = 2;
end

fileNameMaster = [];
fileNameSlave = [];

% Set the input file name
nref = 5;

% generate meshes once forever
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
integration_type = ["SB","EB","RBF"];  % SB, RBF, EB
nInt = 5;
nGP = 3;

rbf_type = 'gauss';



nGrids = 5;
brokenL2 = zeros(nGrids,1);
brokenH1 = zeros(nGrids,1);
h = zeros(nGrids,1);

for i_t = integration_type
  for i = 1:nref
    fprintf('Grid h%i nGP = %i  nInt = %i  Integration: %s \n',i, nGP, nInt, i_t);
    % Import the mesh data into the Mesh object
    fnameBot = fullfile('Mesh','meshes',strcat('bot_',num2str(i),'_',elem_type,'.vtk'));
    fnameTop =  fullfile('Mesh','meshes',strcat('top_',num2str(i),'_',elem_type,'.vtk'));

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

    h(i) = 1/(length(nodesSlave)-1);

    % compute mortar operator
    mortar = Mortar2D(degree,masterMesh,1,slaveMesh,1);
    %D = mortar.D;
    switch i_t
      case 'RBF'
        [D,M] = mortar.computeMortarRBF(6,nInt,rbf_type,mult_type);
      case 'EB'
        [D,M] = mortar.computeMortarElementBased(6,mult_type);
      case 'SB'
        [D,M] = mortar.computeMortarSegmentBased(6,mult_type);
    end

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
      K = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofIs));
        zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofIs)) ;
        KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M';
        zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D';
        zeros(length(dofIs), length(dofM)), zeros(length(dofIs), length(dofS)), -M, D,  zeros(length(dofIs),length(dofIs))];

      f = [f; zeros(length(dofIs),1)];
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

    %     u_anal = @(x,y) 16*x.*y.*(1-x).*(1-y);
%     dux =  @(x,y) 16*y - 16*y.^2 - 32*x.*y + 32*x.*y.^2;
%     duy = @(x,y) 16*x - 16*x.^2 - 32*x.*y + 32*y.*x.^2;

    % computing error
    [L2m,H1m] = computeError(masterMesh,elemsMaster,u_master,u_anal,dux,duy);
    [L2s,H1s] = computeError(slaveMesh,elemsSlave,u_slave,u_anal,dux,duy);

    brokenL2(i) = sqrt(L2m^2+L2s^2);
    brokenH1(i) = sqrt(H1m^2+H1s^2);

    if i == 2 && strcmp(i_t,'RBF') && strcmp(elem_type,'hexa')
      % absolute error of solution
      c = masterMesh.coordinates;
      u_anal_master = u_anal(c(:,1),c(:,2));
      err_rel_master = (u_master-u_anal_master)./u_master;
      err_rel_master(u_master==0) = 0;
      plotFunction(masterMesh,'OUT_errMaster',err_rel_master);
      c = slaveMesh.coordinates;
      u_anal_slave = u_anal(c(:,1),c(:,2));
      err_rel_slave = (u_slave - u_anal_slave)./u_slave;
      err_rel_slave(u_slave==0) = 0;
      plotFunction(slaveMesh,'OUT_errSlave',err_rel_slave);
    end

  end % end refinement loop

  outName = i_t+"_"+elem_type;
  if strcmp(i_t,'RBF')
    outName = outName + "_" + rbf_type;
  end
  out.L2 = brokenL2;
  out.H1 = brokenH1;
  save(outName+".mat",'-struct',"out");
end % end integration type loop







function [L2err,H1err] = computeError(mesh,element,u_app,f_anal,dfx,dfy)
L2err = 0;
H1err = 0;
for el = 1:mesh.nSurfaces
  vtkId = mesh.surfaceVTKType(el);
  elem = getElement(element,vtkId);
  N = getBasisFinGPoints(elem);
  nodes = mesh.surfaces(el,:);
  nodeCoord = mesh.coordinates(nodes,1:2);
  [gradN,dJW] = getDerBasisFAndDet(elem,nodeCoord);
  c_gp = getGPointsLocation(elem,el);
  u_gp = N*u_app(nodes);
  err = u_gp - f_anal(c_gp(:,1),c_gp(:,2));
  err_2 = err.^2;        % squared value of error
  L2errLoc = sum(err_2.*reshape(dJW,[],1));

  grad_uh = pagemtimes(gradN,u_app(nodes));
  grad_uh = squeeze(permute(grad_uh,[3 1 2]));

  grad_u = [dfx(c_gp(:,1),c_gp(:,2)),...
            dfy(c_gp(:,1),c_gp(:,2))];

  grad_err = grad_uh - grad_u;
  grad_err2 = sum(grad_err.^2,2);
  semiH1errLoc = sum(grad_err2.*reshape(dJW,[],1));

  L2err = L2err + L2errLoc;
  H1err = H1err + L2errLoc + semiH1errLoc;
end

L2err = sqrt(L2err);
H1err = sqrt(H1err);
end

