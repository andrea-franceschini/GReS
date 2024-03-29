clear; close all;

% Testing 2D mortar method with RL-RBF
% Mechanical linear elastic problem with 2 non conforming domains and possibly non
% conforming intermediate interface;

% STEP 1: Definition of the geometry. Two separate pieces need to be read
% from GMSH. Need to know the tolopogy of the two interfaces to perform
% immediately the mortar operator computation;

% STEP 2: computing matrix blocks of mortar discretization (use a simple
% function based on the code already written)

% STEP 3: compute the remaining blocks with standard assembly procedure

% STEP 4: output? I should create a paraview file using of the two
% different meshes? or simply a matlab colormap? ask Andrea the easiest
% choice

% DEFINE MODEL
model = ModelType("Poromechanics_FEM");
% Just because it is input for some class constructors

% IMPORT MESHES
topMesh = Mesh();
bottomMesh = Mesh();

% Elastic properties and constitutive tensor
E = 100000;
nu = 0.3;
D = zeros(3);
D([1 5]) = 1-nu;
D([2 4]) = nu;
D(9) = 0.5*(1-2*nu);
D = (E/((1+nu)*(1-2*nu)))*D;
Ku = (E*(1-nu))/((1+nu)*(1-2*nu));
% external force
F = -10000;
% gauss class for 1D element integration
nGP = 4;
gauss = Gauss(12,nGP,1);

fileNameTop = [];
fileNameBottom = [];

% Set the input file name
for k = 1:5
    fileNameTop = [fileNameTop; strcat('Mesh_flat_convergenceTest/TopBlock_tetra_h',num2str(k),'.msh')];
    fileNameBottom = [fileNameBottom; strcat('Mesh_flat_convergenceTest/BottomBlock_tetra_h',num2str(k),'.msh')];
end

flagTop = 'slave';
h = zeros(5,1);
ini_str = strcat("2D MORTAR METHOD: top coarser mesh as ",flagTop, '\n');
fprintf(ini_str)

int_str = ["SB", "RBF", "EB"]; % RBF or SB
for integration = int_str
    if strcmp(integration,'RBF')
        fopen('L2_rbf.txt','w');
        fopen('H1_rbf.txt','w');
    elseif strcmp(integration, 'SB')
        fopen('L2_sb.txt','w');
        fopen('H1_sb.txt','w');
    elseif strcmp(integration, 'EB')
        fopen('L2_eb.txt','w');
        fopen('H1_eb.txt','w');
    end

    for mCount = 1:4
        p_str = strcat(integration,' integration - Mesh size h', num2str(mCount), ' \n');
        fprintf(p_str)
        % Import the mesh data into the Mesh object
        topMesh.importGMSHmesh(fileNameTop(mCount,:));
        bottomMesh.importGMSHmesh(fileNameBottom(mCount,:));

        % assign domains to mortar or slave tag
        % decide if the top domain is master or slave

        if strcmp(flagTop, 'master')
            bottom = 'slave';
            masterMesh = topMesh;
            slaveMesh = bottomMesh;
        elseif strcmp(flagTop, 'slave')
            bottom = 'master';
            slaveMesh = topMesh;
            masterMesh = bottomMesh;
        end

        % Element class for further stiffness matrix computation
        elemsMaster = Elements(masterMesh);
        elemsSlave = Elements(slaveMesh);

        % computing stiffness matrix on each domain
        % DOMAIN 1
        KMaster = stiff(masterMesh, elemsMaster, D);
        KSlave = stiff(slaveMesh, elemsSlave, D);
        % get id of nodes belonging to master and slave interfaces
        nodesMaster = unique(masterMesh.edges(masterMesh.edgeTag == 1,:));
        nodesSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 1,:));

        % compute mortar operator
        if strcmp(integration,'RBF')
            Etmp = compute_mortar(masterMesh, slaveMesh, gauss, 1, 1);
        elseif strcmp(integration, 'SB')
            Etmp = compute_mortar_SB(masterMesh, slaveMesh, 1, 1, 3);
        elseif strcmp(integration, 'EB')
            Etmp = compute_mortar_EB(masterMesh, slaveMesh, gauss, 1, 1);
        end

        % reordering the matrix of the system
        %
        % | dofM = nodi interni master      |
        % | dofS = nodi interni slave       |
        % | dofIm = nodi interfaccia master |

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
        E = zeros(2*nS,2*nM);
        E(2*(1:nS)'-1,2*(1:nM)'-1) = Etmp;
        E(2*(1:nS)',2*(1:nM)') = Etmp;

        K = [Kmm, zeros(length(dofM),length(dofS)), KmIm;
            zeros(length(dofS),length(dofM)), Kss, KsIs*E;
            KmIm', E'*KsIs', KImIm+E'*KIsIs*E];

        % creating global dof array and initializing forcing vector
        listDofs = [dofM';dofS';dofIm'];
        f = zeros(length(listDofs),1);

        % ------------------- APPLY BCS -------------------------------
        % apply boundary conditions to each portion (penalty method)
        nrows = size(K,1);
        penParm = 1e11*max(K,[],'all');

        %------------------- TOP LOAD BCS -----------------------------
        % get Loaded dofs on top edge
        nodesLoad = unique(topMesh.edges(topMesh.edgeTag == 2,:));
        h(mCount) = (length(nodesLoad)-1)^-1;
        % special treatment of extreme points (having force 1/2)
        n_ext = nodesLoad(ismember(topMesh.coordinates(nodesLoad,1),[0; 1]));
        loadDoFext = 2*n_ext;
        loadDoFY = 2*nodesLoad(~ismember(nodesLoad, n_ext));
        loadDoFext = getGlobalDofs(loadDoFext, dofM,dofS,dofIm, flagTop);
        loadDoFY = getGlobalDofs(loadDoFY, dofM,dofS,dofIm, flagTop);
        f(loadDoFext) = 0.5*F/(length(nodesLoad)-1);
        f(loadDoFY) = F/(length(nodesLoad)-1);

        %------------------- BOTTOM FIXED BCS -----------------------------
        % get fixed dofs on bottom edge
        % get nodes
        dirNod = unique(bottomMesh.edges(bottomMesh.edgeTag == 2,:));
        % get associated dofs (in x and y directions)
        dirBotDoF = getDoF(dirNod');
        dirBotDoF = getGlobalDofs(dirBotDoF,dofM,dofS,dofIm,bottom);
        [K,f] = applyDir(dirBotDoF, zeros(length(dirBotDoF),1), K, f);


        %------------------- LATERAL CONSTRAINT  --------------------------

        % LATERAL EDGES: note, only master nodes lying on the interface must be
        % fixed!
        % get nodes from master domain (not in the interface)
        nodesLatMaster = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
        % extract nodes not belonging to the interface
        nodesLatMaster = nodesLatMaster(~ismember(nodesLatMaster, nodesMaster));
        % get corresponding DoFs in the linear system
        dofLatMaster = 2*nodesLatMaster-1; % x direction is fixed
        dofLatMaster = getGlobalDofs(dofLatMaster,dofM,dofS,dofIm,'master');
        % Apply penalty method
        [K,f] = applyDir(dofLatMaster, zeros(length(dofLatMaster),1), K, f);


        % get nodes in the master interface
        nodesLatMaster = unique(masterMesh.edges(masterMesh.edgeTag == 3,:));
        nodesLatInt = nodesLatMaster(ismember(nodesLatMaster, nodesMaster));
        dofLatInt = 2*nodesLatInt-1;
        dofLatInt = getGlobalDofs(dofLatInt,dofM,dofS,dofIm,'interface');
        % Apply penalty method
        [K,f] = applyDir(dofLatInt, zeros(length(dofLatInt),1), K, f);


        % get nodes in the slave domain (not in the interface)
        % get nodes from master domain (not in the interface)
        nodesLatSlave = unique(slaveMesh.edges(slaveMesh.edgeTag == 3,:));
        % extract nodes not belonging to the interface
        nodesLatSlave = nodesLatSlave(~ismember(nodesLatSlave, nodesSlave));
        % get corresponding DoFs in the linear system
        dofLatSlave = 2*nodesLatSlave-1; % x direction is fixed
        dofLatSlave = getGlobalDofs(dofLatSlave,dofM,dofS,dofIm,'slave');
        % Apply penalty method
        [K,f] = applyDir(dofLatSlave, zeros(length(dofLatSlave),1), K, f);

        % ---------------------- SOLVE SYSTEM --------------------------


        % solve linear system
        u = K\f;
        u_slave = E*u(length(dofM)+length(dofS)+1:end);

        % plot solution (use the plot_function already used for the RBF stand alone
        % tests!)

        u_top = zeros(topMesh.nNodes,1);
        u_bottom = zeros(bottomMesh.nNodes,1);

        %collect displacement of master domain and slave domain, according to user assignment;
        if strcmp(flagTop, 'master')
            u_top(dofM) = u(1:length(dofM));
            u_top(dofIm) = u(length(dofM)+length(dofS)+1:end);
            u_bottom(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
            u_bottom(dofIs) = u_slave;
        elseif strcmp(flagTop, 'slave')
            u_bottom(dofM) = u(1:length(dofM));
            u_bottom(dofIm) = u(length(dofM)+length(dofS)+1:end);
            u_top(dofS) = u(length(dofM)+1:length(dofM)+length(dofS));
            u_top(dofIs) = u_slave;
        end

        % Plot results to Paraview

        % fNameTop = strcat('disp_top_h',num2str(mCount),'_',integration);
        % fNameBottom = strcat('disp_bottom_h',num2str(mCount),'_',integration);
        % plotParaview(topMesh,fNameTop, u_top', 'all')
        % plotParaview(bottomMesh,fNameBottom, u_bottom', 'all')


        % Analytical displacements
        u_anal_top = zeros(2*topMesh.nNodes,1);
        u_anal_bot = zeros(2*bottomMesh.nNodes,1);
        u_anal_top(2:2:end) = (F/Ku)*topMesh.coordinates(:,2);
        u_anal_bot(2:2:end) = (F/Ku)*bottomMesh.coordinates(:,2);

        % % solution error (relative)
        err_top = abs(u_anal_top(2:2:end) - u_top(2:2:end))./abs(u_anal_top(2:2:end));
        err_bottom = abs(u_anal_bot(2:2:end) - u_bottom(2:2:end))./abs(u_anal_bot(2:2:end));
        err_top(isinf(err_top)) = 0;
        err_bottom(isinf(err_bottom)) = 0;
        %
        % fNameTopErr = strcat('err_top_h',num2str(mCount),'_',integration);
        % fNameBottomErr = strcat('err_bottom_h',num2str(mCount),'_',integration);
        % plotParaview(topMesh,fNameTopErr, err_top', 'y')
        % plotParaview(bottomMesh,fNameBottomErr, err_bottom', 'y')

        %% Error analysis

        volNodMaster = zeros(masterMesh.nNodes,1);
        volNodSlave = zeros(slaveMesh.nNodes,1);
        %
        % % L2 error norm
        %
        % master mesh
        for el = 1:masterMesh.nSurfaces
            top = masterMesh.surfaces(el, 1:masterMesh.surfaceNumVerts(el));
            volNodMaster(top) = volNodMaster(top) + elemsMaster.tri.findVolume(el)/masterMesh.surfaceNumVerts(el);
        end
        %
        % slave mesh
        for el = 1:slaveMesh.nSurfaces
            top = slaveMesh.surfaces(el, 1:slaveMesh.surfaceNumVerts(el));
            volNodSlave(top) = volNodSlave(top) + elemsSlave.tri.findVolume(el)/slaveMesh.surfaceNumVerts(el);
        end

        if strcmp(flagTop, 'master')
            L2Master = ((u_anal_top(2:2:end) - (u_top(2:2:end)))./u_anal_top(2:2:end)).^2;
            L2Slave = ((u_anal_bot(2:2:end) - (u_bottom(2:2:end)))./u_anal_bot(2:2:end)).^2;
        else
            L2Slave = ((u_anal_top(2:2:end) - (u_top(2:2:end)))./u_anal_top(2:2:end)).^2;
            L2Master = ((u_anal_bot(2:2:end) - (u_bottom(2:2:end)))./u_anal_bot(2:2:end)).^2;
        end
        L2Slave(isinf(L2Slave)) = 0;
        L2Master(isinf(L2Master)) = 0;
        L2Slave(isnan(L2Slave)) = 0;
        L2Master(isnan(L2Master)) = 0;
        L2Slave = sqrt(L2Slave'*volNodSlave);
        L2Master = sqrt(L2Master'*volNodMaster);
        brokenL2 = sqrt(L2Master^2 + L2Slave^2);

        % ERROR IN ENERGY NORM
        % for each element, compute (epsilon-epsilon_h)'*D*(epsilon-epsilon_h)
        % epsilon is the strain tensor computed on each cell (cell-wise
        % computation)
        H1master = zeros(masterMesh.nSurfaces,1);
        H1slave = zeros(slaveMesh.nSurfaces,1);

        % Master surface
        for el = 1:masterMesh.nSurfaces
            top = masterMesh.surfaces(el, 1:masterMesh.surfaceNumVerts(el));
            N = getDerBasisF(elemsMaster.tri,el);
            vol = findVolume(elemsMaster.tri,el);
            B = zeros(3,6);
            B(elemsMaster.indB2D(1:12,2)) = N(elemsMaster.indB2D(1:12,1));
            if strcmp(flagTop,'master')
                e = B*u_anal_top(getDoF(top));
                e_h = B*u_top(getDoF(top));
            else
                e = B*u_anal_bot(getDoF(top));
                e_h = B*u_bottom(getDoF(top));
            end
            H1master(el) = (e-e_h)'*D*(e-e_h)*vol;
        end

        % Slave surface
        for el = 1:slaveMesh.nSurfaces
            top = slaveMesh.surfaces(el, 1:slaveMesh.surfaceNumVerts(el));
            N = getDerBasisF(elemsSlave.tri,el);
            vol = findVolume(elemsSlave.tri,el);
            B = zeros(3,6);
            B(elemsSlave.indB2D(1:12,2)) = N(elemsSlave.indB2D(1:12,1));
            if strcmp(flagTop,'master')
                e = B*u_anal_bot(getDoF(top));
                e_h = B*u_bottom(getDoF(top));
            else
                e = B*u_anal_top(getDoF(top));
                e_h = B*u_top(getDoF(top));
            end
            H1slave(el) = (e-e_h)'*D*(e-e_h)*vol;
        end

        brokenH1 = sqrt(sum(H1master)+sum(H1slave));
        H1slave = sqrt(sum(H1slave));
        H1master = sqrt(sum(H1master));

        % SAVING OUTPUT DATAS IN TEXT FILES
        if strcmp(integration,'RBF')
            fID1 = fopen('L2_rbf.txt', 'a');
            fID2 = fopen('H1_rbf.txt', 'a');
        elseif strcmp(integration, 'SB')
            fID1 = fopen('L2_sb.txt', 'a');
            fID2 = fopen('H1_sb.txt', 'a');
        elseif strcmp(integration, 'EB')
            fID1 = fopen('L2_eb.txt', 'a');
            fID2 = fopen('H1_eb.txt', 'a');
        end
        fprintf(fID1,'%2.5f %2.10f %2.10f %2.10f \n', h(mCount), brokenL2, L2Master, L2Slave);
        fprintf(fID2,'%2.5f %2.10f %2.10f %2.10f \n', h(mCount), brokenH1, H1master, H1slave);
    end
end

%% plot convergence rates

fprintf('Plotting graphs \n')
% H1 error plot
leg_str = [];
figure(1)
if any(strcmp(int_str, 'RBF'))
H1_RBF = load('H1_rbf.txt');
loglog(H1_RBF(:,1), H1_RBF(:,2), '-ro')
hold on
% loglog(H1_RBF(:,1), H1_RBF(:,3), '--r*')
% loglog(H1_RBF(:,1), H1_RBF(:,4), '--rs')
leg_str = [leg_str, "BN RBF integration"];
end


if any(strcmp(int_str, 'SB'))
H1_SB = load('H1_sb.txt');
loglog(H1_SB(:,1), H1_SB(:,2), '-ko')
hold on
% loglog(H1_SB(:,1), H1_SB(:,3), '--k*')
% loglog(H1_SB(:,1), H1_SB(:,4), '--ks')
leg_str = [leg_str, "BN SB integration"];
end

if any(strcmp(int_str, 'EB'))
H1_EB = load('H1_eb.txt');
loglog(H1_EB(:,1), H1_EB(:,2), '-bo')
hold on
% loglog(H1_EB(:,1), H1_EB(:,3), '--b*')
% loglog(H1_EB(:,1), H1_EB(:,4), '--bs')
leg_str = [leg_str, "BN EB integration"];
end
legend(leg_str);
tit_str = strcat("H1 error plot with coarser top domain as ",flagTop);
title(tit_str);


% L2 error plot
leg_str = [];
figure(2)
if any(strcmp(int_str, 'RBF'))
L2_RBF = load('L2_rbf.txt');
loglog(L2_RBF(:,1), L2_RBF(:,2), '-ro')
hold on
%loglog(L2_RBF(:,1), L2_RBF(:,3), '--r*')
%loglog(L2_RBF(:,1), L2_RBF(:,4), '--rs')
leg_str = [leg_str, "BN RBF integration"];
end


if any(strcmp(int_str, 'SB'))
L2_SB = load('L2_sb.txt');
loglog(L2_SB(:,1), L2_SB(:,2), '-ko')
hold on
%loglog(L2_SB(:,1), L2_SB(:,3), '--k*')
%loglog(L2_SB(:,1), L2_SB(:,4), '--ks')
leg_str = [leg_str, "BN SB integration"];
end

if any(strcmp(int_str, 'EB'))
L2_EB = load('L2_eb.txt');
loglog(L2_EB(:,1), L2_EB(:,2), '-bo')
hold on
%loglog(L2_EB(:,1), L2_EB(:,3), '--b*')
%loglog(L2_EB(:,1), L2_EB(:,4), '--bs')
leg_str = [leg_str, "BN EB integration"];
end
legend(leg_str);
tit_str = strcat("L2 error plot with coarser top domain as ",flagTop);
title(tit_str);
fprintf('COMPLETED SUCCESSFULLY \n')