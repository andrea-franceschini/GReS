clc
clear
close all
% Mortar solution of a simple Poisson problem on a unit square domain
% exact solution is: u_ex = 16xy(1-x)(1-y);
% the rhs is f = 32x(1-x) + 32y(1-y);

% flagPlot: true --> results to Paraview in this directory
fPlot = true;

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

degree = 1;
% file Names for input meshes

% selecting master and slave domain
flagTop = 'master';
if strcmp(flagTop, 'master')
    fileNameMaster = 'Mesh/TopBlock_curve.msh';
    fileNameSlave = 'Mesh/BottomBlock_curve.msh';
    strTop = 'masterTop';
elseif strcmp(flagTop, 'slave')
    for k = 1:nGrids
        fileNameSlave = 'Mesh/TopBlock_curve.msh';
        fileNameMaster = 'Mesh/BottomBlock_curve.msh';
    end
    strTop = 'slaveTop';
end

% selecting integration approach
integration = 'EB';  % SB, RBF, EB
nInt = 4;
nGP = 4;
tmp = strcat(flagTop,'TOP');


% Import the mesh data into the Mesh object
masterMesh.importGMSHmesh(fileNameMaster);
slaveMesh.importGMSHmesh(fileNameSlave);
% Element class for further stiffness matrix computation
elemsMaster = Elements(masterMesh);
elemsSlave = Elements(slaveMesh);

% computing Stiffness matrix on top and bottom domainall
[KMaster, aNmaster] = stiffPoisson(masterMesh, elemsMaster);
[KSlave, aNslave] = stiffPoisson(slaveMesh, elemsSlave);

% get id of nodes belonging to master and slave interfaces
nodesMaster = unique(masterMesh.edges(masterMesh.edgeTag == 1,:));
nodesSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 1,:));


% get ID of interface nodes belonging to Dirichlet boundary
% Constant basis functions are considered for slave elements containing
% these nodes (actually this doesn't seem to affect the results)
% boundInt = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
% boundInt = boundInt(ismember(boundInt, nodesSlave));

% compute mortar operator
% and matrices
[E, M, D] = compute_mortar(masterMesh, slaveMesh, [], nInt, nGP, 1, 1, integration, degree);


% reordering the matrix of the system
%
% | dofM = nodi interni master      |
% | dofS = nodi interni slave       |
% | dofIm = nodi interfaccia master |
% | dofIs = nodi interfaccia slave  |
% dof Lagrange multiplier not needed

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
fAn = @(x,y) 32*(x.*(1-x)+y.*(1-y));
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
        zeros(length(dofIs), length(dofM)), zeros(length(dofIs), length(dofS)), -M, D,  zeros(length(dofIs), length(dofIs))];

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
u_anal = @(x,y) 16*x.*y.*(1-x).*(1-y);
% analytical solution
u_anal_master = u_anal(masterMesh.coordinates(:,1), masterMesh.coordinates(:,2));
u_anal_slave = u_anal(slaveMesh.coordinates(:,1), slaveMesh.coordinates(:,2));

postProcMaster = postProc(masterMesh, u_master, u_anal_master);
postProcSlave = postProc(slaveMesh, u_slave, u_anal_slave);

% % Plotting point-wise error
err_rel_master = computeRelError(postProcMaster);
err_rel_slave = computeRelError(postProcSlave);
%
if fPlot
    fNameMaster = strcat(strTop,'_SolMaster','_',integration);
    fNameSlave = strcat(strTop,'_SolSlave','_',integration);
    plotParaview(masterMesh,fNameMaster, u_master', 'x')
    plotParaview(slaveMesh,fNameSlave, u_slave', 'x')
    fNameMaster = strcat(strTop,'_errMaster','_',integration);
    fNameSlave = strcat(strTop,'_errSlave','_',integration);
    plotParaview(masterMesh,fNameMaster, err_rel_master', 'x')
    plotParaview(slaveMesh,fNameSlave, err_rel_slave', 'x')
end

%% Error analysis
L2Master = computeL2error(postProcMaster);
L2Slave = computeL2error(postProcSlave);
brokenL2 = sqrt(L2Master^2 + L2Slave^2);

% H1 error
H1Master = computeH1error(postProcMaster);
H1Slave = computeH1error(postProcSlave);
% Master surface
brokenH1= sqrt(H1Master^2 + H1Slave^2);



% create output structure (or append to existing one) and write to file
% if ~isfile("Results.mat")
%     outStruct= struct('IntegrationType', integration, 'nGP', nGP, 'nInt', nInt,...
%         'MeshSizes', h, 'L2', brokenL2, 'H1', brokenH1);
%     save("Results.mat", "outStruct");
% else
%     out = load('Results.mat', "outStruct");
%     outStruct = out.outStruct;
%     newStruct= struct('IntegrationType', integration, 'nGP', nGP, 'nInt', nInt,...
%         'MeshSizes', h, 'L2', brokenL2, 'H1', brokenH1);
%     outStruct = [outStruct; newStruct];
%     save("Results.mat", "outStruct");
% end

%% PLOT CONVERGENCE PROFILES using struct datas
% out = load('Results.mat', "outStruct");
% outStruct = out.outStruct;
% outStruct = outStruct(1:end);
% results = struct2cell(outStruct);
% 
% 
% tiledlayout(2,1)
% axL2 = nexttile;
% axH1 = nexttile;
% axL2.NextPlot = "add";
% axH1.NextPlot = "add";
% 
% 
% param = 'nInt';
% marker = {'o', '^', 's', '*'};
% getFirst = @(v)v{1}; 
% getprop = @(options, idx)getFirst(circshift(options,-idx+1));
% 
% % plot convergence profiles
% for ii = 1:size(results,2)
%     int = results{1,ii};
%     if strcmp(param, 'GP')
%         parm = results{2,ii};
%     else strcmp(param,'nInt')
%         parm = results{3,ii};
%     end
%     switch int
%         case 'RBF'
%             lgd = strcat(int,' - ',param, ' = ',num2str(parm));
%             color = 'b';
%             mark = getprop(marker,ii);
%         case 'SB'
%             lgd = strcat(int);
%             color = 'r';
%             mark = 'o';
%         case 'EB'
%             lgd = strcat(int);
%             color = 'g';
%             mark = 'o';
%     end
%     loglog(axL2,results{4,ii}, results{5,ii},'Color',color,'Marker',...
%         mark,'DisplayName',lgd,'LineWidth', 1);
%     loglog(axH1,results{4,ii}, results{6,ii},'Color',color,'Marker',...
%         mark,'DisplayName',lgd,'LineWidth',1);
% end
% set(axL2, 'YScale', 'log')
% set(axL2, 'XScale', 'log')
% set(axH1, 'YScale', 'log')
% set(axH1, 'XScale', 'log')
% legend(axL2);
% legend(axH1);
% xlabel('Mesh size');
% ylabel(axL2,'L2 norm of error')
% ylabel(axH1,'H1 norm of error')
% 
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 12);
% a = get(axL2,'XTickLabel');
% set(axL2,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 9)
% a = get(axH1,'XTickLabel');
% set(axH1,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 9)