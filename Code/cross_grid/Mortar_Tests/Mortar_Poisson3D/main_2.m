clear
close all
% Mortar solution of a simple Poisson problem on a 3D domain
% exact solution is: u_ex = cos(pi*y)*cos(pi*z)*(2-2x+pi*cos(pi*x));
% the rhs is f = ; cos(pi*y)*cos(pi*z)*(-2 -3*pi^2*sin(pi*x) - 4*pi^2*x+ 2*pi^2*x^2
%

% Analytical solution
u_anal = @(x,y,z) cos(pi*y).*cos(pi*z).*(2*x - x.^2 + sin(pi*x));

% Analtical rhs
h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
fAn = @(x,y,z) -cos(pi*y).*cos(pi*z).*h(x);

% selecting solution method
% COND --> condansated approach
% SP --> saddle point matrix
sol_scheme = 'COND';

% Parameters
nGrids = 4;
fPlot = true;
nGP = 2;
nInt = 4;
type = 'gauss';
brokenL2 = zeros(nGrids,1);
brokenH1 = zeros(nGrids,1);
h = zeros(nGrids,1);

% IMPORT MESHES
masterMesh = Mesh();
slaveMesh = Mesh();

% leftMesh.importGMSHmesh('mesh/domainLeft_H1.msh');
% rightMesh.importGMSHmesh('mesh/domainRight_H1.msh');

% file Names for input meshes

fileNameMaster = [];
fileNameSlave = [];

% fileNameMaster = 'meshCurve/LeftBlock_curve.msh';
% fileNameSlave = 'meshCurve/RightBlock_curve.msh';


% Set the input file name
% selecting master and slave domain
flagLeft = 'slave';
if strcmp(flagLeft, 'master')
    for k = 1:nGrids
        fileNameMaster = [fileNameMaster; strcat('meshConv/LeftBlock_h',num2str(k),'.msh')];
        fileNameSlave = [fileNameSlave; strcat('meshConv/RightBlock_h',num2str(k),'.msh')];
    end
    strLeft = 'masterLeft';
elseif strcmp(flagLeft, 'slave')
    for k = 1:nGrids
        fileNameSlave = [fileNameSlave; strcat('MeshConv/LeftBlock_h',num2str(k),'.msh')];
        fileNameMaster = [fileNameMaster; strcat('MeshConv/RightBlock_h',num2str(k),'.msh')];
    end
    strLeft = 'slaveLeft';
end

% Gauss integration for stiffness matrix computation
gauss = Gauss(12,2,3);

for mCount = 1:nGrids
    fprintf('Grid h%i nGP = %i  nInt = %i \n',mCount, nGP, nInt);
    % fprintf(p_str)
    % % Import the mesh data into the Mesh object
    masterMesh.importGMSHmesh(fileNameMaster(mCount,:));
    slaveMesh.importGMSHmesh(fileNameSlave(mCount,:));
    % Element class for further stiffness matrix computation
    elemsMaster = Elements(masterMesh, gauss);
    elemsSlave = Elements(slaveMesh, gauss);

    % computing Stiffness matrix on top and bottom domain
    [KMaster, vNmaster] = stiffPoisson3D(masterMesh, elemsMaster);
    [KSlave, vNslave] = stiffPoisson3D(slaveMesh, elemsSlave);

    mortar = Mortar3D(1,masterMesh,1,slaveMesh,1);
    % compute mortar operator and matrices
    %[D,M,~,tRBF,E] = mortar.computeMortarRBF(nGP,nInt,type);
    switch type
        case 'gauss'
            [D,M,tRBF,~,E] = mortar.computeMortarRBF(nGP,nInt,type);
        case 'eb'
            [D,M,tEB,~,E] = mortar.computeMortarElementBased(nGP);
    end
    % reordering the matrix of the system
    %
    % | dofM = nodi interni master      |
    % | dofS = nodi interni slave       |
    % | dofIm = nodi interfaccia master |
    % | dofIs = nodi interfaccia slave  |
    % dof Lagrange multiplier not needed

    dofIm = mortar.nodesMaster;
    dofM = (1:masterMesh.nNodes)';
    dofM = dofM(~ismember(dofM,dofIm));
    dofIs = mortar.nodesSlave;
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
    zM = [masterMesh.coordinates(dofM,3);
        masterMesh.coordinates(dofIm,3)];
    xS = [slaveMesh.coordinates(dofS,1);
        slaveMesh.coordinates(dofIs,1)];
    yS = [slaveMesh.coordinates(dofS,2);
        slaveMesh.coordinates(dofIs,2)];
    zS = [slaveMesh.coordinates(dofS,3);
        slaveMesh.coordinates(dofIs,3)];
    % compute forcing term (LM rows are set to 0)
    x = [xM;xS];
    y = [yM;yS];
    z = [zM;zS];
    f = fAn(x,y,z);
    % compute forcing vector and areanod
    f = f.*[vNmaster(dofM);vNmaster(dofIm);vNslave(dofS);vNslave(dofIs)];
    % reorder system rhs
    f = [f(1:length(dofM));
        f(length(dofM)+length(dofIm)+1:length(dofM)+length(dofIm)+length(dofS));
        f(length(dofM)+1:length(dofM)+length(dofIm));
        f(length(dofM)+length(dofIm)+length(dofS)+1:end)];

    if strcmp(sol_scheme, 'SP')
        % complete saddle point matrix
        K = [Kmm, sparse(length(dofM),length(dofS)), KmIm, sparse(length(dofM),length(dofIs)), sparse(length(dofM),length(dofIs));
            sparse(length(dofS),length(dofM)), Kss, sparse(length(dofS),length(dofIm)), KsIs, sparse(length(dofS),length(dofIs));
            KmIm', sparse(length(dofIm),length(dofS)), KImIm, sparse(length(dofIm),length(dofIs)), -M';
            sparse(length(dofIs),length(dofM)), KsIs', sparse(length(dofIs),length(dofIm)), KIsIs, D';
            sparse(length(dofIs),length(dofM)), sparse(length(dofIs),length(dofS)), -M, D, sparse(length(dofIs),length(dofIs))];
        %
        f = [f; zeros(length(dofIs),1)];
        listDofs = [dofM;dofS;dofIm; dofIs; dofIs];
    else
        K = [Kmm, sparse(length(dofM),length(dofS)), KmIm;
            sparse(length(dofS),length(dofM)), Kss, KsIs*E;
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
    nodesLatMaster = unique(masterMesh.surfaces(masterMesh.surfaceTag == 2,:));
    % extract nodes not belonging to the interface
    nodesLatMaster = nodesLatMaster(~ismember(nodesLatMaster, mortar.nodesMaster));
    % get corresponding DoFs in the linear system
    dofLatMaster = nodesLatMaster;
    dofLatMaster = getGlobalDofs(dofLatMaster,dofM,dofS,dofIm,dofIs,'master');
    % evalutate analytical solution on dirichlet nodes of the master domain
    x = masterMesh.coordinates(nodesLatMaster,1);
    y = masterMesh.coordinates(nodesLatMaster,2);
    z = masterMesh.coordinates(nodesLatMaster,3);
    fDir = u_anal(x,y,z);
    
    % Apply Dirichlet BCs with standard method (1 in the main diagonal)
    [K,f] = applyDir(dofLatMaster, fDir, K, f);

    % slave domain
    nodesLatSlave = unique(slaveMesh.surfaces(slaveMesh.surfaceTag == 2,:));
    % extract nodes not belonging to the interface
    nodesLatSlave = nodesLatSlave(~ismember(nodesLatSlave, mortar.nodesSlave));
    % get corresponding DoFs in the linear system
    dofLatSlave = nodesLatSlave; % x direction is fixed
    dofLatSlave = getGlobalDofs(dofLatSlave,dofM,dofS,dofIm,dofIs,'slave');
    % evaluate analytical solution on dirichlet nodes of the slave domain
    x = slaveMesh.coordinates(nodesLatSlave,1);
    y = slaveMesh.coordinates(nodesLatSlave,2);
    z = slaveMesh.coordinates(nodesLatSlave,3);
    fDir = u_anal(x,y,z);

    % Apply Dirichlet BCs
    [K,f] = applyDir(dofLatSlave, fDir, K, f);

    % master interface
    nodesLatMaster = unique(masterMesh.surfaces(masterMesh.surfaceTag == 2,:));
    nodesLatInt = nodesLatMaster(ismember(nodesLatMaster, mortar.nodesMaster));
    dofLatInt = nodesLatInt;
    dofLatInt = getGlobalDofs(dofLatInt,dofM,dofS,dofIm,dofIs,'interfaceMaster');
    % evalutate analytical solution on dirichlet nodes of the master
    % interface
    x = masterMesh.coordinates(nodesLatInt,1);
    y = masterMesh.coordinates(nodesLatInt,2);
    z = masterMesh.coordinates(nodesLatInt,3);
    fDir = u_anal(x,y,z);
    % Apply Dirichlet BCs
    [K,f] = applyDir(dofLatInt, fDir, K, f);
    %
    % if strcmp(sol_scheme,'SP')
    %     % slave interface
    %     nodesLatSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
    %     nodesLatInt = nodesLatSlave(ismember(nodesLatSlave, mortar.nodesSlave));
    %     dofLatInt = nodesLatInt;
    %     dofLatInt = getGlobalDofs(dofLatInt,dofM,dofS,dofIm,dofIs,'interfaceSlave');
    %     % Apply Dirichlet BCs
    %     [K,f] = applyDir(dofLatInt, zeros(length(dofLatInt),1), K, f);
    % end


    % -------------------------------- SOLVE SYSTEM --------------------------

    % solve linear system
    u = K\f;
    u_master = zeros(masterMesh.nNodes,1);
    u_slave = zeros(slaveMesh.nNodes,1);
    l = length(dofM)+length(dofS)+length(dofIm);
    if strcmp(sol_scheme,'SP')
        u_s = u(l+1:l+length(dofIs));
    else
        u_s = E*u(length(dofM)+length(dofS)+1:l);
    end

    %collect displacement of master domain and slave domain
    u_master(dofM) = u(1:length(dofM));
    u_master(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
    u_slave(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
    u_slave(dofIs) = u_s;

    u_anal_master = u_anal(masterMesh.coordinates(:,1), masterMesh.coordinates(:,2), masterMesh.coordinates(:,3));
    u_anal_slave = u_anal(slaveMesh.coordinates(:,1), slaveMesh.coordinates(:,2),  slaveMesh.coordinates(:,3));

    % Plotting point-wise error
    postProcMaster = postProc(masterMesh, u_master, u_anal_master,gauss);
    postProcSlave = postProc(slaveMesh, u_slave, u_anal_slave,gauss);

    % % Plotting point-wise error
    err_rel_master = computeRelError(postProcMaster);
    err_rel_slave = computeRelError(postProcSlave);

    if fPlot ==  true
        fNameMaster = strcat(strLeft,'_SolMaster','_h',num2str(mCount));
        fNameSlave = strcat(strLeft,'_SolSlave','_h',num2str(mCount));
        plotParaview(masterMesh,fNameMaster, u_master', 'x')
        plotParaview(slaveMesh,fNameSlave, u_slave', 'x')

        % fNameMaster = strcat(strLeft,'_errMaster','_h',num2str(mCount));
        % fNameSlave = strcat(strLeft,'_errSlave','_h',num2str(mCount));
        % plotParaview(masterMesh,fNameMaster, err_rel_master', 'x')
        % plotParaview(slaveMesh,fNameSlave, err_rel_slave', 'x')
    end

    % if strcmp(fPlot,'curve')
    %     fNameMaster = strcat(strLeft,'_SolMasterCurve','_h',num2str(mCount));
    %     fNameSlave = strcat(strLeft,'_SolSlaveCurve','_h',num2str(mCount));
    %     plotParaview(masterMesh,fNameMaster, u_master', 'x')
    %     plotParaview(slaveMesh,fNameSlave, u_slave', 'x')
    % 
    %     fNameMaster = strcat(strLeft,'_errMasterCurve','_h',num2str(mCount));
    %     fNameSlave = strcat(strLeft,'_errSlaveCurve','_h',num2str(mCount));
    %     plotParaview(masterMesh,fNameMaster, err_rel_master', 'x')
    %     plotParaview(slaveMesh,fNameSlave, err_rel_slave', 'x')
    % end


    % L2 error
    L2Master = computeL2error(postProcMaster);
    L2Slave = computeL2error(postProcSlave);
    brokenL2(mCount) = sqrt(L2Master^2 + L2Slave^2);

    % H1 error
    H1Master = computeH1error(postProcMaster);
    H1Slave = computeH1error(postProcSlave);
    % Master surface
    brokenH1(mCount) = sqrt(H1Master^2 + H1Slave^2);
    h(mCount) = getGridSize(postProcMaster);
end

%% Save results in text file
switch type
    case 'gauss'
        name = strcat(flagLeft,'L2_',type,'_Int',num2str(nInt));
        fID = fopen(strcat('Results\',name,'.dat'),'w');
        fprintf(fID,'%2.6e \n',brokenL2);

        name = strcat(flagLeft,'H1_',type,'_Int',num2str(nInt));
        fID = fopen(strcat('Results\',name,'.dat'),'w');
        fprintf(fID,'%2.6e \n',brokenH1);
    case 'eb'
        name = strcat(flagLeft,'L2_eb');
        fID = fopen(strcat('Results\',name,'.dat'),'w');
        fprintf(fID,'%2.6e \n',brokenL2);

        name = strcat(flagLeft,'H1_eb');
        fID = fopen(strcat('Results\',name,'.dat'),'w');
        fprintf(fID,'%2.6e \n',brokenH1);
end

   


% create output structure (or appen to existing one) and write to file
% if ~isfile("Results.mat")
%     outStruct= struct('nGP', nGP, 'nInt', nInt,...
%         'MeshSizes', h, 'L2', brokenL2, 'H1', brokenH1);
%     save("Results.mat", "outStruct");
% else
%     out = load('Results.mat', "outStruct");
%     outStruct = out.outStruct; 
%     newStruct= struct('nGP', nGP, 'nInt', nInt,...
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
% tiledlayout(2,1)
% axL2 = nexttile;
% axH1 = nexttile;
% axL2.NextPlot = "add";
% axH1.NextPlot = "add";
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
%         parm = results{1,ii};
%     else strcmp(param,'nInt')
%         parm = results{2,ii};
%     end
%     lgd = strcat(param,' = ',num2str(parm));
%     color = 'k';
%     mark = getprop(marker,ii);
%     loglog(axL2,results{3,ii}, results{4,ii},'Color',color,'Marker',mark,...
%         'DisplayName',lgd,'LineWidth', 1);
%     loglog(axH1,results{3,ii}, results{5,ii},'Color',color,'Marker',mark,...
%         'DisplayName',lgd,'LineWidth', 1);
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