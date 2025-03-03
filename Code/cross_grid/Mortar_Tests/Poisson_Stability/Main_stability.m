clc
clear
close all
% Mortar solution of a simple Poisson problem on a unit square domain
% exact solution is: u_ex = 16xy(1-x)(1-y);
% the rhs is f = 32x(1-x) + 32y(1-y);

% flagPlot: true --> results to Paraview in this directory
fPlot = true;


% IMPORT MESHES
masterMesh = Mesh();
slaveMesh = Mesh();

% Import meshes (no convergence profiles)
% topMesh.importGMSHmesh('Mesh/TopBlock_curve.msh');
% bottomMesh.importGMSHmesh('Mesh/BottomBlock_curve.msh');

degree = 1;
% file Names for input meshes

fileNameMaster = [];
fileNameSlave = [];
mult_type = 'dual';

% Set the input file name
nGrids = 3;
% selecting master and slave domain

% selecting integration approach
for type = ["DUAL","CONF","P0"]
   nInt = 4;
   nGP = 4;
   switch type
      case 'DUAL'
         brokenL2_RBF = zeros(nGrids,1);
         brokenH1_RBF = zeros(nGrids,1);
         leg = 'Dual';
      case 'P0'
         brokenL2_P0 = zeros(nGrids,1);
         brokenH1_P0 = zeros(nGrids,1);
         leg = 'P0 stabilized';
      case 'CONF'
         brokenL2_CONF = zeros(nGrids,1);
         brokenH1_CONF = zeros(nGrids,1);
         leg = 'Conforming (FINE)';
      case 'UNBIASED'
         leg = 'Unbiased (Puso stabilization)';
   end
   h = zeros(nGrids,1);
   %grids = 2; % row vector for selecting in which to compute pointwise error
   NX0 = 4; % Initial number of elements in the top mesh
   rat = 2;
   for mCount = 1:nGrids
      NelX = NX0*(2^(mCount-1)); 
      fprintf('Grid h%i nGP = %i  nInt = %i  Integration: %s \n',mCount, nGP, nInt, type);
      % Import the mesh data into the Mesh object
      getMeshPois2D('Mesh_conv/domain_top.geo','top',NelX);
      if strcmp(type,'CONF')
         getMeshPois2D('Mesh_conv/domain_top.geo','top',rat*NelX);
      end
      getMeshPois2D('Mesh_conv/domain_bottom.geo','bottom',rat*NelX);
      masterMesh = Mesh(); slaveMesh = Mesh();
      masterMesh.importGMSHmesh('Mesh_conv/top.msh');
      slaveMesh.importGMSHmesh('Mesh_conv/bottom.msh');
      % Element class for further stiffness matrix computation
      elemsMaster = Elements(masterMesh);
      elemsSlave = Elements(slaveMesh);
      % instance of 2D mortar class
      mortar = Mortar2D(1,masterMesh,1,slaveMesh,1);
      h(mCount) = 1/NelX;
      %D = mortar.D;
      switch type
         case 'DUAL'
            [D,M] = mortar.computeMortarRBF(nGP,nInt,'gauss',mult_type);
         case 'P0'
            [D,M] = mortar.computeMortarConstant(nGP,nInt);
         case 'CONF'
            [D,M] = mortar.computeConfCrossGridMat();
         case 'UNBIASED'
            mortar1 = Mortar2D(1,masterMesh,1,slaveMesh,1);
            mortar2 = Mortar2D(1,slaveMesh,1,masterMesh,1);
            [D1,M1] = mortar1.computeMortarRBF(nGP,4,'gauss',mult_type);
            [D2,M2] = mortar2.computeMortarRBF(nGP,4,'gauss',mult_type);
      end
      dof = DofMap(masterMesh,mortar.nodesMaster,slaveMesh,mortar.nodesSlave);
      dofIm = DofMap.getCompDoF(mortar.nodesMaster,1);
      dofM = DofMap.getCompDoF((1:masterMesh.nNodes)',1);
      dofM = dofM(~ismember(dofM,dofIm));
      dofIs = DofMap.getCompDoF(mortar.nodesSlave,1);
      dofS = DofMap.getCompDoF((1:slaveMesh.nNodes)',1);
      dofS = dofS(~ismember(dofS,dofIs));
      switch type
         case 'P0'
            dofMult = DofMap.getCompDoF((1:mortar.nElSlave)',1);
         case {'DUAL','CONF'}
            E = D\M;
            % handle end points
            % D(3,:) = D(3,:) + D(1,:);
            % D(end,:) = D(end,:) + D(2,:);
            % M(3,:) = M(3,:) + M(1,:);
            % M(end,:) = M(end,:) + M(2,:);
            % D(1,:) = []; D(1,:) = [];
            % M(1,:) = []; M(1,:) = [];
            dofMult = DofMap.getCompDoF(mortar.nodesSlave(1:end),1);
      end

      % computing Stiffness matrix on top and bottom domainall
      [KMaster, aNmaster] = stiffPoisson(masterMesh, elemsMaster);
      [KSlave, aNslave] = stiffPoisson(slaveMesh, elemsSlave);
      % get id of nodes belonging to master and slave interfaces
      Kmm = KMaster(dofM,dofM);
      KmIm = KMaster(dofM,dofIm);
      Kss = KSlave(dofS,dofS);
      KsIs = KSlave(dofS, dofIs);
      KImIm = KMaster(dofIm,dofIm);
      KIsIs = KSlave(dofIs,dofIs);

      hM = h(mCount);
      hS = 1/mortar.nElSlave;
      if strcmp(type,'P0')
         H = hM*mortar.computePressureJumpMat();
      else
         H = zeros(length(dofMult),length(dofMult));
      end

      %H = zeros(length(dofMult),length(dofMult));

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
      %fAn = @(x,y) 32*(x.*(1-x)+y.*(1-y));
      fAn = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
      f = fAn(x,y);
      % compute forcing vector and areanod
      f = f.*[aNmaster(dofM);aNmaster(dofIm);aNslave(dofS);aNslave(dofIs)];
      % reorder system rhs
      f = [f(1:length(dofM));
         f(length(dofM)+length(dofIm)+1:length(dofM)+length(dofIm)+length(dofS));
         f(length(dofM)+1:length(dofM)+length(dofIm));
         f(length(dofM)+length(dofIm)+length(dofS)+1:end)];

      % complete saddle point matrix
      switch type
         case 'UNBIASED'
            s = (hM);
            r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofIm)), zeros(length(dofM),length(dofIs))];
            r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofIm)), zeros(length(dofS),length(dofIs))] ;
            r3 =  [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), 0.5*D2', -0.5*M1'];
            r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, -0.5*M2', 0.5*D1'];
            r5 = [zeros(length(dofIm), length(dofM)), zeros(length(dofIm), length(dofS)), 0.5*D2, -0.5*M2,  s*D2, s*M2];
            r6 = [zeros(length(dofIs), length(dofM)), zeros(length(dofIs), length(dofS)), -0.5*M1, 0.5*D1,  s*M1, s*D1];
            K = [r1;r2;r3;r4;r5;r6];
            K = [r1;r2;r3;r4;r5;r6];
            f = [f; zeros(length(dofIm)+length(dofIs),1)];
         otherwise
            r1 = [Kmm, zeros(length(dofM),length(dofS)), KmIm, zeros(length(dofM),length(dofIs)), zeros(length(dofM),length(dofMult))];
            r2 = [zeros(length(dofS),length(dofM)), Kss, zeros(length(dofS),length(dofIm)), KsIs, zeros(length(dofS),length(dofMult))] ;
            r3 = [KmIm', zeros(length(dofIm),length(dofS)), KImIm, zeros(length(dofIm),length(dofIs)), -M'];
            r4 = [zeros(length(dofIs), length(dofM)), KsIs', zeros(length(dofIs), length(dofIm)), KIsIs, D'];
            r5 = [zeros(length(dofMult), length(dofM)), zeros(length(dofMult), length(dofS)), -M, D, -H];
            K = [r1;r2;r3;r4;r5];
            f = [f; zeros(length(dofMult),1)];
      end

      % ------------------------------ APPLY BCS -------------------------------
      % homogeneous Dirichlet BCs
      % master domain
      % get nodes from master domain (not in the interface)
      nodesLatMaster = unique(masterMesh.edges(masterMesh.edgeTag == 2,:));
      % extract nodes not belonging to the interface
      %nodesLatMaster = nodesLatMaster(~ismember(nodesLatMaster, mortar.nodesMaster));
      % get corresponding DoFs in the linear system
      dofLatMaster = dof.getDoF(nodesLatMaster,'master',1);
      % Apply Dirichlet BC
      [K,f] = applyDir(dofLatMaster, zeros(length(dofLatMaster),1), K, f);
      % slave domain
      nodesLatSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
      %nodesLatSlave = nodesLatSlave(~ismember(nodesLatSlave, mortar.nodesSlave));
      %
      % get corresponding DoFs in the linear system
      dofLatSlave = dof.getDoF(nodesLatSlave,'slave',1);
      % Apply Dirichlet BCs
      [K,f] = applyDir(dofLatSlave, zeros(length(dofLatSlave),1), K, f);

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
      u_s = u(l+1:l+length(dofIs));

      u_slave(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
      u_slave(dofIs) = u_s;


      %% Plotting solutions
      u_anal = @(x,y) sin(pi*x).*sin(pi*y);
      % analytical solution
      u_anal_master = u_anal(masterMesh.coordinates(:,1), masterMesh.coordinates(:,2));
      u_anal_slave = u_anal(slaveMesh.coordinates(:,1), slaveMesh.coordinates(:,2));

      postProcMaster = postProc(masterMesh, u_master, u_anal_master);
      postProcSlave = postProc(slaveMesh, u_slave, u_anal_slave);

      %% Error analysis
      L2Master = computeL2error(postProcMaster);
      L2Slave = computeL2error(postProcSlave);
      % H1 error
      H1Master = computeH1error(postProcMaster);
      H1Slave = computeH1error(postProcSlave);
      %
      bL2 = sqrt(L2Master^2 + L2Slave^2);
      bH1 = sqrt(H1Master^2 + H1Slave^2);
      switch type
         case 'DUAL'
            brokenL2_RBF(mCount) = bL2;
            brokenH1_RBF(mCount) = bH1;
         case 'P0'
            brokenL2_P0(mCount) = bL2;
            brokenH1_P0(mCount) = bH1;
         case 'CONF'
            brokenL2_CONF(mCount) = bL2;
            brokenH1_CONF(mCount) = bH1;
      end
   end

   mult = u(end-numel(dofMult)+1:end);


   % plot multipliers along the interface

   switch type
      case {'DUAL','CONF','UNBIASED'}
         [x,id] = sort(slaveMesh.coordinates(mortar.nodesSlave(3:end),1));
      case 'P0'
         c = [slaveMesh.coordinates(mortar.slaveTopol(:,1),1), slaveMesh.coordinates(mortar.slaveTopol(:,2),1)];
         x = 0.5*sum(c,2);
         [x,id] = sort(x);
   end
   figure(1)
   plot(x,mult(id),'s-','LineWidth',1,'DisplayName',leg)
   xlabel('x-coordinate interface')
   ylabel('Multiplier')
   hold on
      %Stabilization trough projection
   if strcmp(type,'DUAL')
      % expand the solution to deleted multipliers
      mult2 = zeros(numel(dofMult)+2,1);
      mult2(1) = mult(1);
      mult2(2) = mult(end);
      mult2(3:end) = mult;
      mortarSwap = Mortar2D(1,slaveMesh,1,masterMesh,1);
      [Dsw, Msw] = mortarSwap.computeMortarSegmentBased(2,mult_type);
      Esw = Dsw\Msw;
      mult = E*(Esw*mult2);
      [x,id] = sort(slaveMesh.coordinates(mortar.nodesSlave,1));
      plot(x,mult(id),'s-','LineWidth',1,'DisplayName','Naive projection')
   end

   if fPlot
      plotParaview(masterMesh,strcat('SOLmaster_',type,'_h',num2str(mCount)), u_master', 'x')
      plotParaview(slaveMesh,strcat('SOLslave_',type,'_h',num2str(mCount)), u_slave', 'x')
      plotParaview(masterMesh,strcat('ERRmaster_',type,'_h',num2str(mCount)), abs(u_anal_master'-u_master'), 'x')
      plotParaview(slaveMesh,strcat('ERRslave_',type,'_h',num2str(mCount)), abs(u_anal_slave'-u_slave'), 'x')
   end
end

figure(1)
% plot([0 1],[0 0],'k-','LineWidth',1,'Exact')
legend show

figure(2)
tiledlayout(1,1) % Crea un layout 1x2

% Primo plot: Broken L2 norm
nexttile
loglog(h, brokenL2_RBF, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'DUAL - Broken L2 Norm') % RBF in blu
hold on
loglog(h, brokenL2_P0, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'P0 - Broken L2 Norm') % P0 in rosso
loglog(h, brokenL2_CONF, 'k-d', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'CONFORMING - Broken L2 Norm') % CONF in nero

% Linea di riferimento h^2 spostata verso il basso
loglog(h, 0.5 * h.^2 * (brokenL2_RBF(1)/h(1)^2), 'K-.', 'LineWidth', 1, 'DisplayName', 'h^2')

% Miglioramenti grafici
legend('Location', 'southwest', 'FontSize', 12) % Sposta la legenda
xlabel('h', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Error', 'FontSize', 12, 'FontWeight', 'bold')
title('Broken L2 Norm', 'FontSize', 14, 'FontWeight', 'bold')

% Secondo plot: Broken H1 norm
% nexttile
% loglog(h, brokenH1_RBF, 'b-o', 'DisplayName', 'DUAL - Broken H1 Norm') % RBF in blu
% hold on
% loglog(h, brokenH1_P0, 'r-s', 'DisplayName', 'P0 - Broken H1 Norm') % P0 in rosso
% loglog(h, brokenH1_CONF, 'k-s', 'DisplayName', 'CONFORMING - Broken L2 Norm') % P0 in rosso
% loglog(h, h * (brokenH1_RBF(1)/h(1)), 'k-.', 'DisplayName', 'h') % Riferimento
% legend show
% xlabel('h')
% ylabel('Error')
% title('Broken H1 Norm')
