clear
close all
% Mortar solution of a simple Poisson problem on a unit square domain
% exact solution is: u_ex = 16xy(1-x)(1-y);
% the rhs is f = 32x(1-x) + 32y(1-y);

% flagPlot: true --> results to Paraview in this directory
fPlot = false;

% selecting solution method
% COND --> condansated approach
% SP --> saddle point matrix
sol_scheme = 'SP';

% IMPORT MESHES
topMesh = Mesh();
bottomMesh = Mesh();

topMesh.importGMSHmesh('Mesh/TopBlock_tetra.msh');
bottomMesh.importGMSHmesh('Mesh/BottomBlock_tetra.msh');

% file Names for input meshes

fileNameTop = [];
fileNameBottom = [];

% Set the input file name
for k = 1:4
    fileNameTop = [fileNameTop; strcat('Mesh_conv/TopBlock_tetra_h',num2str(k),'.msh')];
    fileNameBottom = [fileNameBottom; strcat('Mesh_conv/BottomBlock_tetra_h',num2str(k),'.msh')];
end

% selecting master and slave domain
flagTop = 'master';
if strcmp(flagTop, 'master')
    bottom = 'slave';
    masterMesh = topMesh;
    slaveMesh = bottomMesh;
elseif strcmp(flagTop, 'slave')
    bottom = 'master';
    slaveMesh = topMesh;
    masterMesh = bottomMesh;
end

% selecting integration approach
int_str = ["SB","RBF","EB"];

tmp = strcat(flagTop,'TOP');
% open files for writing convergence results
for integration = int_str
    if strcmp(integration,'RBF')
        fopen(strcat(tmp,'L2_rbf.txt'),'w');
        fopen(strcat(tmp,'H1_rbf.txt'),'w');
    elseif strcmp(integration, 'SB')
        fopen(strcat(tmp,'L2_sb.txt'),'w');
        fopen(strcat(tmp,'H1_sb.txt'),'w');
    elseif strcmp(integration, 'EB')
        fopen(strcat(tmp,'L2_eb.txt'),'w');
        fopen(strcat(tmp,'H1_eb.txt'),'w');
    end


    nGrids = 4;

    for mCount = 1:nGrids
        p_str = strcat(integration,' integration - Mesh size h', num2str(mCount), ' \n');
        fprintf(p_str)
        % Import the mesh data into the Mesh object
        topMesh.importGMSHmesh(fileNameTop(mCount,:));
        bottomMesh.importGMSHmesh(fileNameBottom(mCount,:));
        % Element class for further stiffness matrix computation
        elemsMaster = Elements(masterMesh);
        elemsSlave = Elements(slaveMesh);

        % computing Stiffness matrix on top and bottom domainall
        [KMaster, aNmaster] = stiffPoisson(masterMesh, elemsMaster);
        [KSlave, aNslave] = stiffPoisson(slaveMesh, elemsSlave);

        % get id of nodes belonging to master and slave interfaces
        nodesMaster = unique(masterMesh.edges(masterMesh.edgeTag == 1,:));
        nodesSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 1,:));

        h = 1/(length(nodesSlave)-1);

        % get ID of interface nodes belonging to Dirichlet boundary
        % Constant basis functions are considered for slave elements containing
        % these nodes (actually this doesn't seem to affect the results)
        boundInt = unique(slaveMesh.edges(slaveMesh.edgeTag == 2,:));
        boundInt = boundInt(ismember(boundInt, nodesSlave));

        % compute mortar operator and matrices
        if strcmp(integration,'RBF')
            [E, M, D] = compute_mortar(masterMesh, slaveMesh, [], 20, 6, 1, 1, "RBF");
        elseif strcmp(integration, 'SB')
            [E, M, D] = compute_mortar(masterMesh, slaveMesh, [], 13, 3, 1, 1, "SB");
        elseif strcmp(integration, 'EB')
            [E, M, D] = compute_mortar(masterMesh, slaveMesh, [], 13, 6, 1, 1, "EB");
            %[Etmp] = compute_mortar_EB(masterMesh, slaveMesh, gaussEB, 1, 1);
        end

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

        u_top = zeros(topMesh.nNodes,1);
        u_bottom = zeros(bottomMesh.nNodes,1);

        %collect displacement of master domain and slave domain, according to user assignment;
        if strcmp(flagTop, 'master')
            u_top(dofM) = u(1:length(dofM));
            u_top(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
            u_bottom(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
            u_bottom(dofIs) = u_slave;
        elseif strcmp(flagTop, 'slave')
            u_bottom(dofM) = u(1:length(dofM));
            u_bottom(dofIm) = u(length(dofM)+length(dofS)+1:length(dofM)+length(dofS)+length(dofIm));
            u_top(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
            u_top(dofIs) = u_slave;
        end

        % Error analysis
        u_anal = @(x,y) 16*x.*y.*(1-x).*(1-y);
        % analytical solution
        u_anal_top = u_anal(topMesh.coordinates(:,1), topMesh.coordinates(:,2));
        u_anal_bot = u_anal(bottomMesh.coordinates(:,1), bottomMesh.coordinates(:,2));

        % Plotting point-wise error
        err_top_rel = abs((u_anal_top - u_top)./u_anal_top);
        err_bot_rel = abs((u_anal_bot - u_bottom)./u_anal_bot);
        err_top_rel(isinf(err_top_rel)) = 0;
        err_bot_rel(isinf(err_bot_rel)) = 0;
        err_top_rel(isnan(err_top_rel)) = 0;
        err_bot_rel(isnan(err_bot_rel)) = 0;
        %
        err_top = abs(u_anal_top - u_top);
        err_bot = abs(u_anal_bot - u_bottom);

        if fPlot
            % Plot results to Paraview
            if strcmp(flagTop,'master')
                strTop = 'masterTop';
            else
                strTop = 'slaveTop';
            end

            fNameTop = strcat(strTop,'_SolTop','_',integration,'_h',num2str(mCount));
            fNameBottom = strcat(strTop,'_SolBot','_',integration,'_h',num2str(mCount));
            plotParaview(topMesh,fNameTop, u_top', 'x')
            plotParaview(bottomMesh,fNameBottom, u_bottom', 'x')

            fNameTop = strcat(strTop,'_errTop','_',integration,'_h',num2str(mCount));
            fNameBottom = strcat(strTop,'_errBot','_',integration,'_h',num2str(mCount));
            plotParaview(topMesh,fNameTop, err_top_rel', 'x')
            plotParaview(bottomMesh,fNameBottom, err_bot_rel', 'x')
        end

        % Error norms
        % L2 error
        if strcmp(flagTop, 'master')
            L2Master = (u_anal_top - u_top).^2;
            L2Slave = (u_anal_bot - u_bottom).^2;
        else
            L2Slave = (u_anal_top - u_top).^2;
            L2Master = (u_anal_bot - u_bottom).^2;
        end

        L2Slave(isinf(L2Slave)) = 0;
        L2Master(isinf(L2Master)) = 0;
        L2Slave(isnan(L2Slave)) = 0;
        L2Master(isnan(L2Master)) = 0;
        L2Slave = sqrt(L2Slave'*aNslave);
        L2Master = sqrt(L2Master'*aNmaster);
        brokenL2 = sqrt(L2Master^2 + L2Slave^2);

        % H1 error
        H1Master = zeros(masterMesh.nSurfaces,1);
        H1Slave = zeros(slaveMesh.nSurfaces,1);
        % Master surface
        for el = 1:masterMesh.nSurfaces
            top = masterMesh.surfaces(el, 1:masterMesh.surfaceNumVerts(el));
            N = getDerBasisF(elemsMaster.tri,el);
            vol = findVolume(elemsMaster.tri,el);
            if strcmp(flagTop,'master')
                u_ex = u_anal_top(top);
                u_h = u_top(top);
            else
                u_ex = u_anal_bot(top);
                u_h = u_bottom(top);
            end
            H1Master(el) = (N*(u_ex-u_h))'*(N*(u_ex-u_h));
            H1Master(el) = H1Master(el)*vol;
        end

        % Slave surface
        for el = 1:slaveMesh.nSurfaces
            top = slaveMesh.surfaces(el, 1:slaveMesh.surfaceNumVerts(el));
            N = getDerBasisF(elemsSlave.tri,el);
            vol = findVolume(elemsSlave.tri,el);
            if strcmp(flagTop,'master')
                u_ex = u_anal_bot(top);
                u_h = u_bottom(top);
            else
                u_ex = u_anal_top(top);
                u_h = u_top(top);
            end
            H1Slave(el) = (N*(u_ex-u_h))'*(N*(u_ex-u_h));
            H1Slave(el) = H1Slave(el)*vol;
        end
        
        %H1 seminorms squared
        H1Slave = sqrt(sum(H1Slave));
        H1Master = sqrt(sum(H1Master));
        % H1Slave = sqrt(H1Slave+L2Slave^2);
        % H1Master = sqrt(H1Master+L2Master^2);
        brokenH1 = sqrt(H1Master^2 + H1Slave^2);

        % SAVING OUTPUT DATAS IN TEXT FILES
        if strcmp(integration,'RBF')
            fID1 = fopen(strcat(tmp,'L2_rbf.txt'), 'a');
            fID2 = fopen(strcat(tmp,'H1_rbf.txt'), 'a');
        elseif strcmp(integration, 'SB')
            fID1 = fopen(strcat(tmp,'L2_sb.txt'), 'a');
            fID2 = fopen(strcat(tmp,'H1_sb.txt'), 'a');
        elseif strcmp(integration, 'EB')
            fID1 = fopen(strcat(tmp,'L2_eb.txt'), 'a');
            fID2 = fopen(strcat(tmp,'H1_eb.txt'), 'a');
        end
        fprintf(fID1,'%2.5f %2.10f %2.10f %2.10f \n', h, brokenL2, L2Master, L2Slave);
        fprintf(fID2,'%2.5f %2.10f %2.10f %2.10f \n', h, brokenH1, H1Master, H1Slave);
    end
end

%% PLOT CONVERGENCE PROFILES

fprintf('Plotting graphs \n')
% H1 error plot
leg_str = [];
figure(1)
if any(strcmp(int_str, 'RBF'))
H1_RBF = load(strcat(tmp,'H1_rbf.txt'));
loglog(H1_RBF(:,1), H1_RBF(:,2), '--ro', 'LineWidth', 1, 'MarkerSize', 8.5)
hold on
%loglog(H1_RBF(:,1), H1RBFgp4, '--r^', 'LineWidth', 1, 'MarkerSize', 8.5)
% loglog(H1_RBF(:,1), H1RBFgp6, '--rs',  'LineWidth', 1, 'MarkerSize', 8.5)
leg_str = [leg_str, "RBF integration"];
end


if any(strcmp(int_str, 'SB'))
H1_SB = load(strcat(tmp,'H1_sb.txt'));
loglog(H1_SB(:,1), H1_SB(:,2), '-ko',  'LineWidth', 1, 'MarkerSize', 8.5)
hold on
% loglog(H1_SB(:,1), H1_SB(:,3), '--k*')
% loglog(H1_SB(:,1), H1_SB(:,4), '--ks')
leg_str = [leg_str, "SB integration"];
end

if any(strcmp(int_str, 'EB'))
H1_EB = load(strcat(tmp,'H1_eb.txt'));
loglog(H1_EB(:,1), H1_EB(:,2), '-bo', 'LineWidth', 1, 'MarkerSize', 8.5)
hold on
% loglog(H1_EB(:,1), H1_EB(:,3), '--b*')
% loglog(H1_EB(:,1), H1_EB(:,4), '--bs')
leg_str = [leg_str, "EB integration"];
end
legend(leg_str, 'Location', 'northwest');
tit_str = strcat("H1 error plot with coarser top domain as ",flagTop);
title(tit_str);
xticks([])
xlabel('Mesh size')
ylabel('H1 energy broken error seminorm')
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 12)

% L2 error plot
leg_str = [];
figure(2)
if any(strcmp(int_str, 'RBF'))
L2_RBF = load(strcat(tmp,'L2_rbf.txt'));
loglog(L2_RBF(:,1), L2_RBF(:,2), '--ro', 'LineWidth', 1, 'MarkerSize', 8.5)
hold on
%loglog(H1_RBF(:,1), H1RBFgp4, '--r^', 'LineWidth', 1, 'MarkerSize', 8.5)
% loglog(H1_RBF(:,1), H1RBFgp6, '--rs',  'LineWidth', 1, 'MarkerSize', 8.5)
leg_str = [leg_str, "RBF integration"];
end


if any(strcmp(int_str, 'SB'))
L2_SB = load(strcat(tmp,'L2_sb.txt'));
loglog(L2_SB(:,1), L2_SB(:,2), '-ko',  'LineWidth', 1, 'MarkerSize', 8.5)
hold on
% loglog(H1_SB(:,1), H1_SB(:,3), '--k*')
% loglog(H1_SB(:,1), H1_SB(:,4), '--ks')
leg_str = [leg_str, "SB integration"];
end

if any(strcmp(int_str, 'EB'))
L2_EB = load(strcat(tmp,'L2_eb.txt'));
loglog(L2_EB(:,1), L2_EB(:,2), '-bo', 'LineWidth', 1, 'MarkerSize', 8.5)
hold on
% loglog(H1_EB(:,1), H1_EB(:,3), '--b*')
% loglog(H1_EB(:,1), H1_EB(:,4), '--bs')
leg_str = [leg_str, "EB integration"];
end
legend(leg_str, 'Location', 'northwest');
tit_str = strcat("L2 error plot with coarser top domain as ",flagTop);
title(tit_str);
xticks([])
xlabel('Mesh size')
ylabel('L2 broken error norm')
set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 14);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 12)















