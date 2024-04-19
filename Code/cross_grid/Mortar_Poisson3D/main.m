clear
close all
% Mortar solution of a simple Poisson problem on a 3D domain
% exact solution is: u_ex = cos(pi*y)*cos(pi*z)*(2-2x+pi*cos(pi*x));
% the rhs is f = ; cos(pi*y)*cos(pi*z)*(-2 -3*pi^2*sin(pi*x) - 4*pi^2*x+ 2*pi^2*x^2
%

% flagPlot: true --> results to Paraview in this directory
fPlot = true;

% selecting solution method
% COND --> condansated approach
% SP --> saddle point matrix
sol_scheme = 'COND';

% IMPORT MESHES
leftMesh = Mesh();
rightMesh = Mesh();

leftMesh.importGMSHmesh('mesh/domainLeft_H1.msh');
rightMesh.importGMSHmesh('mesh/domainRight_H1.msh');

% file Names for input meshes

% fileNameTop = [];
% fileNameBottom = [];


% selecting master and slave domain
flagLeft = 'master';
if strcmp(flagLeft, 'master')
    right = 'slave';
    masterMesh = leftMesh;
    slaveMesh = rightMesh;
elseif strcmp(flagLeft, 'slave')
    bottom = 'master';
    masterMesh = rightMesh;
    slaveMesh = leftMesh;
end

nGrids = 1;

% Gauss integration for stiffness matrix computation
gauss = Gauss(12,3,3);

for mCount = 1:nGrids
    % p_str = strcat(integration,' integration - Mesh size h', num2str(mCount), ' \n');
    % fprintf(p_str)
    % % Import the mesh data into the Mesh object
    % topMesh.importGMSHmesh(fileNameTop(mCount,:));
    % bottomMesh.importGMSHmesh(fileNameBottom(mCount,:));
    % Element class for further stiffness matrix computation
    elemsMaster = Elements(masterMesh, gauss);
    elemsSlave = Elements(slaveMesh, gauss);

    % computing Stiffness matrix on top and bottom domainall
    [KMaster, vNmaster] = stiffPoisson3D(masterMesh, elemsMaster);
    [KSlave, vNslave] = stiffPoisson3D(slaveMesh, elemsSlave);

    % get id of nodes belonging to master and slave interfaces
    nodesMaster = unique(masterMesh.surfaces(masterMesh.surfaceTag == 1,:));
    nodesSlave = unique(slaveMesh.surfaces(slaveMesh.surfaceTag == 1,:));

    % get interfaces as mesh objects
    intMaster = masterMesh.getSurfaceMesh(1);
    intSlave = slaveMesh.getSurfaceMesh(1);

    % connectivity matrix between the interfaces
    cs = ContactSearching(intMaster,intSlave,18);

    % h = 1/(length(nodesSlave)-1);

    % get ID of interface nodes belonging to Dirichlet boundary
    % Constant basis functions are considered for slave elements containing
    % these nodes (actually this doesn't seem to affect the results)
    % boundInt = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
    % boundInt = boundInt(ismember(boundInt, nodesSlave));
    % 
    nGPrbf = 5;
    nINTrbf = 10;
    % compute mortar operator and matrices
    [E, M, D] = compute_mortar3D(intMaster, intSlave, cs.elemConnectivity, nGPrbf, nINTrbf);


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
    h = @(x) -2-3*pi^2*sin(pi*x)-4*pi^2*x+2*pi^2*x.^2;
    fAn = @(x,y,z) cos(pi*y).*cos(pi*z).*h(x);
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
    nodesLatMaster = unique(masterMesh.surfaces(masterMesh.surfaceTag == 1,:));
    % extract nodes not belonging to the interface
    nodesLatMaster = nodesLatMaster(~ismember(nodesLatMaster, nodesMaster));
    % get corresponding DoFs in the linear system
    dofLatMaster = nodesLatMaster;
    dofLatMaster = getGlobalDofs(dofLatMaster,dofM,dofS,dofIm,dofIs,'master');
    % Apply Dirichlet BCs with standard method (1 in the main diagonal)
    [K,f] = applyDir(dofLatMaster, zeros(length(dofLatMaster),1), K, f);

    % slave domain
    nodesLatSlave = unique(slaveMesh.surfaces(slaveMesh.surfaceTag == 1,:));
    % extract nodes not belonging to the interface
    nodesLatSlave = nodesLatSlave(~ismember(nodesLatSlave, nodesSlave));
    % get corresponding DoFs in the linear system
    dofLatSlave = nodesLatSlave; % x direction is fixed
    dofLatSlave = getGlobalDofs(dofLatSlave,dofM,dofS,dofIm,dofIs,'slave');
    % Apply Dirichlet BCs
    [K,f] = applyDir(dofLatSlave, zeros(length(dofLatSlave),1), K, f);

    % master interface
    nodesLatMaster = unique(masterMesh.surfaces(masterMesh.surfaceTag == 1,:));
    nodesLatInt = nodesLatMaster(ismember(nodesLatMaster, nodesMaster));
    dofLatInt = nodesLatInt;
    dofLatInt = getGlobalDofs(dofLatInt,dofM,dofS,dofIm,dofIs,'interfaceMaster');
    % Apply Dirichlet BCs
    [K,f] = applyDir(dofLatInt, zeros(length(dofLatInt),1), K, f);
    %
    % if strcmp(sol_scheme,'SP')
    %     % slave interface
    %     nodesLatSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
    %     nodesLatInt = nodesLatSlave(ismember(nodesLatSlave, nodesSlave));
    %     dofLatInt = nodesLatInt;
    %     dofLatInt = getGlobalDofs(dofLatInt,dofM,dofS,dofIm,dofIs,'interfaceSlave');
    %     % Apply Dirichlet BCs
    %     [K,f] = applyDir(dofLatInt, zeros(length(dofLatInt),1), K, f);
    % end




    % -------------------------------- SOLVE SYSTEM --------------------------

    % solve linear system
    u = K\f;

    l = length(dofM)+length(dofS)+length(dofIm);
    if strcmp(sol_scheme,'SP')
        u_slave = u(l+1:l+length(dofIs));
    else
        u_slave = E*u(length(dofM)+length(dofS)+1:l);
    end

    u_left = zeros(leftMesh.nNodes,1);
    u_right = zeros(rightMesh.nNodes,1);

    %collect displacement of master domain and slave domain, according to user assignment;
    if strcmp(flagLeft, 'master')
        u_left(dofM) = u(1:length(dofM));
        u_left(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
        u_right(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
        u_right(dofIs) = u_slave;
    elseif strcmp(flagLeft, 'slave')
        u_right(dofM) = u(1:length(dofM));
        u_right(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
        u_left(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
        u_left(dofIs) = u_slave;
    end




    % Computing errors
    u_anal = @(x,y,z) cos(pi*y).*cos(pi*z).*(2*x - x.^2 + sin(pi*x));
    % analytical solution
    u_anal_left = u_anal(leftMesh.coordinates(:,1), leftMesh.coordinates(:,2), leftMesh.coordinates(:,3));
    u_anal_right = u_anal(rightMesh.coordinates(:,1), rightMesh.coordinates(:,2),  rightMesh.coordinates(:,3));

    % Plotting point-wise error
    err_left_rel = abs((u_anal_left - u_left)./u_anal_left);
    err_right_rel = abs((u_anal_right - u_right)./u_anal_right);
    err_left_rel(isinf(err_left_rel)) = 0;
    err_right_rel(isinf(err_right_rel)) = 0;
    err_left_rel(isnan(err_left_rel)) = 0;
    err_right_rel(isnan(err_right_rel)) = 0;
    %
    err_left = abs(u_anal_left - u_left);
    err_right = abs(u_anal_right - u_right);

    if fPlot
        % Plot results to Paraview
        if strcmp(flagLeft,'master')
            strLeft = 'masterLeft';
        else
            strLeft = 'slaveLeft';
        end

        fNameLeft = strcat(strLeft,'_SolLeft','_h',num2str(mCount));
        fNameRight = strcat(strLeft,'_SolRight','_h',num2str(mCount));
        plotParaview(leftMesh,fNameLeft, u_left', 'x')
        plotParaview(rightMesh,fNameRight, u_right', 'x')


        fNameLeft = strcat(strLeft,'_errLeft','_h',num2str(mCount));
        fNameRight = strcat(strLeft,'_errRight','_h',num2str(mCount));
        plotParaview(leftMesh,fNameLeft, err_left_rel', 'x')
        plotParaview(rightMesh,fNameRight, err_right_rel', 'x')
    end
end

%     %% Error Analysis
%     % L2 error
%     if strcmp(flagTop, 'master')
%         L2Master = (u_anal_top - u_left).^2;
%         L2Slave = (u_anal_bot - u_right).^2;
%     else
%         L2Slave = (u_anal_top - u_left).^2;
%         L2Master = (u_anal_bot - u_right).^2;
%     end
% 
%     L2Slave(isinf(L2Slave)) = 0;
%     L2Master(isinf(L2Master)) = 0;
%     L2Slave(isnan(L2Slave)) = 0;
%     L2Master(isnan(L2Master)) = 0;
%     L2Slave = sqrt(L2Slave'*aNslave);
%     L2Master = sqrt(L2Master'*aNmaster);
%     brokenL2 = sqrt(L2Master^2 + L2Slave^2);
% 
%     % H1 error
% 
%     H1Master = zeros(masterMesh.nSurfaces,1);
%     H1Slave = zeros(slaveMesh.nSurfaces,1);
%     % Master surface
%     for el = 1:masterMesh.nSurfaces
%         top = masterMesh.surfaces(el, 1:masterMesh.surfaceNumVerts(el));
%         N = getDerBasisF(elemsMaster.tri,el);
%         vol = findVolume(elemsMaster.tri,el);
%         if strcmp(flagTop,'master')
%             u_ex = u_anal_top(top);
%             u_h = u_left(top);
%         else
%             u_ex = u_anal_bot(top);
%             u_h = u_right(top);
%         end
%         H1Master(el) = (N*(u_ex-u_h))'*(N*(u_ex-u_h));
%         H1Master(el) = H1Master(el)*vol;
%     end
% 
%     % Slave surface
%     for el = 1:slaveMesh.nSurfaces
%         top = slaveMesh.surfaces(el, 1:slaveMesh.surfaceNumVerts(el));
%         N = getDerBasisF(elemsSlave.tri,el);
%         vol = findVolume(elemsSlave.tri,el);
%         if strcmp(flagTop,'master')
%             u_ex = u_anal_bot(top);
%             u_h = u_right(top);
%         else
%             u_ex = u_anal_top(top);
%             u_h = u_left(top);
%         end
%         H1Slave(el) = (N*(u_ex-u_h))'*(N*(u_ex-u_h));
%         H1Slave(el) = H1Slave(el)*vol;
%     end
% 
%     %H1 seminorms squared
%     H1Slave = sqrt(sum(H1Slave));
%     H1Master = sqrt(sum(H1Master));
%     % H1Slave = sqrt(H1Slave+L2Slave^2);
%     % H1Master = sqrt(H1Master+L2Master^2);
%     brokenH1 = sqrt(H1Master^2 + H1Slave^2);
% 
%     % SAVING OUTPUT DATAS IN TEXT FILES
%     if strcmp(integration,'RBF')
%         fID1 = fopen(strcat(tmp,'L2_rbf.txt'), 'a');
%         fID2 = fopen(strcat(tmp,'H1_rbf.txt'), 'a');
%     elseif strcmp(integration, 'SB')
%         fID1 = fopen(strcat(tmp,'L2_sb.txt'), 'a');
%         fID2 = fopen(strcat(tmp,'H1_sb.txt'), 'a');
%     elseif strcmp(integration, 'EB')
%         fID1 = fopen(strcat(tmp,'L2_eb.txt'), 'a');
%         fID2 = fopen(strcat(tmp,'H1_eb.txt'), 'a');
%     end
%     fprintf(fID1,'%2.5f %2.10f %2.10f %2.10f \n', h, brokenL2, L2Master, L2Slave);
%     fprintf(fID2,'%2.5f %2.10f %2.10f %2.10f \n', h, brokenH1, H1Master, H1Slave);
% end
% 
% %% POST PROCESSING ---- PLOT CONVERGENCE PROFILES
% 
% fprintf('Plotting graphs \n')
% % H1 error plot
% leg_str = [];
% figure(1)
% if any(strcmp(int_str, 'RBF'))
%     H1_RBF = load(strcat(tmp,'H1_rbf.txt'));
%     loglog(H1_RBF(:,1), H1_RBF(:,2), '--ro', 'LineWidth', 1, 'MarkerSize', 8.5)
%     hold on
%     %loglog(H1_RBF(:,1), H1RBFgp4, '--r^', 'LineWidth', 1, 'MarkerSize', 8.5)
%     % loglog(H1_RBF(:,1), H1RBFgp6, '--rs',  'LineWidth', 1, 'MarkerSize', 8.5)
%     leg_str = [leg_str, "RBF integration"];
% end
% 
% 
% if any(strcmp(int_str, 'SB'))
%     H1_SB = load(strcat(tmp,'H1_sb.txt'));
%     loglog(H1_SB(:,1), H1_SB(:,2), '-ko',  'LineWidth', 1, 'MarkerSize', 8.5)
%     hold on
%     % loglog(H1_SB(:,1), H1_SB(:,3), '--k*')
%     % loglog(H1_SB(:,1), H1_SB(:,4), '--ks')
%     leg_str = [leg_str, "SB integration"];
% end
% 
% if any(strcmp(int_str, 'EB'))
%     H1_EB = load(strcat(tmp,'H1_eb.txt'));
%     loglog(H1_EB(:,1), H1_EB(:,2), '-bo', 'LineWidth', 1, 'MarkerSize', 8.5)
%     hold on
%     % loglog(H1_EB(:,1), H1_EB(:,3), '--b*')
%     % loglog(H1_EB(:,1), H1_EB(:,4), '--bs')
%     leg_str = [leg_str, "EB integration"];
% end
% legend(leg_str, 'Location', 'northwest');
% tit_str = strcat("H1 error plot with coarser top domain as ",flagTop);
% title(tit_str);
% xticks([])
% xlabel('Mesh size')
% ylabel('H1 energy broken error seminorm')
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 12)
% 
% % L2 error plot
% leg_str = [];
% figure(2)
% if any(strcmp(int_str, 'RBF'))
%     L2_RBF = load(strcat(tmp,'L2_rbf.txt'));
%     loglog(L2_RBF(:,1), L2_RBF(:,2), '--ro', 'LineWidth', 1, 'MarkerSize', 8.5)
%     hold on
%     %loglog(H1_RBF(:,1), H1RBFgp4, '--r^', 'LineWidth', 1, 'MarkerSize', 8.5)
%     % loglog(H1_RBF(:,1), H1RBFgp6, '--rs',  'LineWidth', 1, 'MarkerSize', 8.5)
%     leg_str = [leg_str, "RBF integration"];
% end
% 
% 
% if any(strcmp(int_str, 'SB'))
%     L2_SB = load(strcat(tmp,'L2_sb.txt'));
%     loglog(L2_SB(:,1), L2_SB(:,2), '-ko',  'LineWidth', 1, 'MarkerSize', 8.5)
%     hold on
%     % loglog(H1_SB(:,1), H1_SB(:,3), '--k*')
%     % loglog(H1_SB(:,1), H1_SB(:,4), '--ks')
%     leg_str = [leg_str, "SB integration"];
% end
% 
% if any(strcmp(int_str, 'EB'))
%     L2_EB = load(strcat(tmp,'L2_eb.txt'));
%     loglog(L2_EB(:,1), L2_EB(:,2), '-bo', 'LineWidth', 1, 'MarkerSize', 8.5)
%     hold on
%     % loglog(H1_EB(:,1), H1_EB(:,3), '--b*')
%     % loglog(H1_EB(:,1), H1_EB(:,4), '--bs')
%     leg_str = [leg_str, "EB integration"];
% end
% legend(leg_str, 'Location', 'northwest');
% tit_str = strcat("L2 error plot with coarser top domain as ",flagTop);
% title(tit_str);
% xticks([])
% xlabel('Mesh size')
% ylabel('L2 broken error norm')
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 12)
