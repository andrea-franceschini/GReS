%%
%TESTING CROSS INTEGRATION OVER A PATCH OF NON MATCHING 1D ELEMENTS
%
% a generic analytical function is defined over the mortar surface
% the code computes the mortar operator E such that:
%               u_s = E * u_m
%                E = inv(D)*M
%
% D is the slave matrix (a typical mass matrix, which is diagonal
% when using a dual lagrange multiplier space)
%
% M is the cross matrix, which is computed with both an "exact" algorithm
% and the RL-RBF procedure
%
% the L2 error norm of the input function and the obtained u_s is computed
% changing the number of interpolation points for the RBF interpolation and
% changing the number of Gauss Points for the element-based discontinous
% integration
% This script uses quadratic basis function for both primal variables and
% Lagrange multipliers
% The nodal based approach with scaling has been dropped 

%% INPUT DATAS

close all; clear

% INPUT
degree = 2;
tol = 1.e-4;
nSizes = 9;
errNormRBF = zeros(nSizes,1);
% master and slave surfaces node position along X axis
nMnodes = 9; % number of nodes for the first mesh refinment
errNormEX = errNormRBF;
errNormEB = errNormRBF;
for sizeCount = 1:nSizes
    nMaster = nMnodes*(2^(sizeCount-1));
    if mod(nMaster,2) == 0
        nMaster = nMaster+1;
    end
    nSlave = round((2/3)*nMaster);
    if mod(nSlave,2) == 0
        nSlave = nSlave+1;
    end
    master = zeros(nMaster,2); % coordinates of master side
    slave = zeros(nSlave,2);
    x1 = -1;
    x2 = 1;
    master(:,1) = (linspace(x1,x2,nMaster))';
    slave(:,1) = (linspace(x1,x2,nSlave))';

    % y-axis (to add possible overlapping bewteen the interfaces)
    % we set the maximum overlapping to 1/10 of the grid size
    sm = abs(master(2)-master(1));
    ss = abs(slave(2)-slave(1));
    fact = 0;
    master(:,2) = -fact*sm*rand(nMaster,1)+fact*sm*rand(nMaster,1);
    slave(:,2) = -fact*ss*rand(nSlave,1)+fact*ss*rand(nSlave,1);


    % number of RBF interpolation points for each element
    nInt = 10;

    % Number of integration points for RBF testing (GP class taken from GReS)
    nGP = 4 ;

    % Inizialize output matrices
    D = zeros(length(slave), length(slave));

    % initializing mortar matrix resulting from RBF computation

    M_RBF = zeros(length(master), length(slave));
    M_EB = zeros(length(master), length(slave));
    M_EX = zeros(size(master,1), size(slave,1));
    % Build a topology matrix for master/slave surfs based on nodes position
    mastertop = build_topol_quad(master(:,1));
    slavetop = build_topol_quad(slave(:,1));

    % Element conntectivity matrix
    elem_connectivity = zeros(size(mastertop,1), size(slavetop,1));
    for i = 1:size(mastertop,1)
        a = master(mastertop(i,1));
        b = master(mastertop(i,3));
        % loop trough master element to find connectivity
        for j = 1:size(slavetop,1)
            c = slave(slavetop(j,1));
            d = slave(slavetop(j,3));
            if ~any([a>d,c>b])
                % intersecting
                elem_connectivity(i,j) = 1;
            end
        end
    end

    %% COMPUTE SLAVE MATRIX D
    % use Gauss integration with 3 Gauss points
    g = Gauss(12,nGP,1);
    gpRef = g.coord;
    gpWeight = g.weight;
    for sID = 1:size(slavetop,1)
        % elem length
        i1 = slavetop(sID,1);
        i2 = slavetop(sID,2);
        i3 = slavetop(sID,3);
        h = norm(slave(i1,:)-slave(i3,:));
        N = computeQuadraticSF(gpRef);
        Dloc = 0.5*h*N'*(N.*gpWeight);
        D([i1 i2 i3], [i1 i2 i3]) = D([i1 i2 i3], [i1 i2 i3]) + Dloc;
    end

    %% Element-wise RBF matrix computation  
    % rbf interpolation is applied to a single element for each of its
    % nodes
    % It's more efficient for different reasons:
    % a smoother function is interpolated (no kinks in the node)
    % smaller linear systems are computed
    % for the rescaling, only one system is solved for each element
    for i = 1:size(mastertop,1)
        % get toy linear shape functions for detecting the intersections
        [valsLin, ptsLin] = computeBasisF1D(i, 4, mastertop, master, degree-2);
        % get actual shape function values and interpolation points
        % coordinates
        [vals, ptsInt] = computeBasisF1D(i, nInt, mastertop, master, degree);
        fiMM = zeros(length(ptsInt), length(ptsInt));
        fiMMLin = zeros(length(ptsLin), length(ptsLin));
        % circumradius of the master element
        r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2);
        % compute interpolation matrix for quadratic function
        for ii = 1:length(ptsInt)
            d = sqrt((ptsInt(:,1) - ptsInt(ii,1)).^2 + ((ptsInt(:,2)) - ptsInt(ii,2)).^2);
            rbf = pos(1-d./r).^4.*(1+d./r);
            fiMM(ii,:) = rbf;
        end
        % compute interpolation matrix for linear functions
        rLin = sqrt((max(ptsLin(:,1)) - min(ptsLin(:,1)))^2 + (max(ptsLin(:,2)) - min(ptsLin(:,2)))^2);
        for ii = 1:length(ptsLin)
            d = sqrt((ptsLin(:,1) - ptsLin(ii,1)).^2 + ((ptsLin(:,2)) - ptsLin(ii,2)).^2);
            rbf = pos(1-d./rLin).^4.*(1+d./rLin);
            fiMMLin(ii,:) = rbf;
        end
        % solve local system to get weight of interpolant
        % compute RBF weight of interpolant
        weightFLin = fiMMLin\valsLin;
        weight1Lin = fiMMLin\ones(size(ptsLin,1),1);
        weightF = fiMM\vals;
        weight1 = fiMM\ones(size(ptsInt,1),1);
        % get slave elements connected to the master
        slave_elems = find(elem_connectivity(i,:));
        % loop trough connected slave elements
        for j = slave_elems
            % get nodes of slave elements
            a = slave(slavetop(j,1),:);
            b = slave(slavetop(j,3),:);
            % get Gauss Points position in the real space
            ptsGauss = ref2nod(gpRef, a, b);
            fiNM = zeros(size(ptsGauss,1), size(ptsInt,1));
            fiNMLin = zeros(size(ptsGauss,1), size(ptsLin,1));
            % compute interpolant of quadratic basis on local Gauss points
            for jj = 1:size(ptsGauss,1)
                d = sqrt((ptsInt(:,1) - ptsGauss(jj,1)).^2 + ((ptsInt(:,2)) - ptsGauss(jj,2)).^2);
                rbf = pos(1-d./r).^4.*(1+d./r);
                fiNM(jj,:) = rbf;
            end
            % compute interpolant of linear basis on local Gauss points
            for jj = 1:size(ptsGauss,1)
                d = sqrt((ptsLin(:,1) - ptsGauss(jj,1)).^2 + ((ptsLin(:,2)) - ptsGauss(jj,2)).^2);
                rbf = pos(1-d./rLin).^4.*(1+d./rLin);
                fiNMLin(jj,:) = rbf;
            end
            Ntmp = (fiNMLin*weightFLin)./(fiNMLin*weight1Lin);
            NMaster = (fiNM*weightF)./(fiNM*weight1);
            %NMaster = fiNM*weightF;
            % automatically detect GP outside the master support using
            % cheap linear interpolation
            id = any([Ntmp < 0, Ntmp > 1, isnan(Ntmp)],2);
            NMaster(id,:) = 0;
            % get basis function on Gauss points
            NSlave = computeQuadraticSF(gpRef);
            % Compute local M matrix using Gauss integration
            Mloc = NMaster'*(NSlave.*gpWeight);
            % Apply the determinant of the jacobian (lenght of the
            % segment)
            Mloc = Mloc*(0.5*sqrt((b(1)-a(1))^2 + (b(2)-a(2))^2));
            % populate cross matrix
            idMaster = [mastertop(i,1); mastertop(i,2); mastertop(i,3)];
            idSlave = [slavetop(j,1); slavetop(j,2); slavetop(j,3)];
            M_RBF(idMaster, idSlave) = M_RBF(idMaster, idSlave) + Mloc;
        end
    end

    %% SEGMENT-BASED CROSS MATRIX COMPUTATION (valid for flat interfaces)
    if all(master(:,2)==0)
        for i = 1:size(mastertop,1)
            % get element of the support
            a = master(mastertop(i,1),:);
            b = master(mastertop(i,3),:);
            % get shape function values and interpolation points coordinate
            slave_elems = find(elem_connectivity(i,:));
            % loop trough connected slave elements
            for j = slave_elems
                % get nodes
                c = slave(slavetop(j,1),:);
                d = slave(slavetop(j,3),:);
                % find intersection between master and slave (just looking at the
                % x-projection)
                tmp = sort([a(1) b(1) c(1) d(1)]);
                i1 = [tmp(2) 0];
                i2 = [tmp(3) 0];
                % get GP in the intersection element
                gPts = ref2nod(gpRef, i1, i2);
                % compute basis function on master element points
                gMaster = nod2ref(gPts, a, b);
                basisMaster = computeQuadraticSF(gMaster);
                % compute basis function on slave elements points
                gSlave = nod2ref(gPts, c, d);
                basisSlave = computeQuadraticSF(gSlave);
                basisSlave = basisSlave.*gpWeight;
                matVal = basisMaster'*basisSlave;
                matVal = matVal*0.5*abs(i2(1)-i1(1));
                idMaster = [mastertop(i,1); mastertop(i,2); mastertop(i,3)];
                idSlave = [slavetop(j,1); slavetop(j,2); slavetop(j,3)];
                M_EX(idMaster, idSlave) = M_EX(idMaster, idSlave) + matVal;
            end
        end
    end

    %% ELEMENT-BASED CROSS MATRIX COMPUTATION (valid for flat interfaces)
    if all(master(:,2)==0)
        gEx = Gauss(12,nGP,1);
        gpRef = gEx.coord;
        gpWeight = gEx.weight;
        for i = 1:size(mastertop,1)
            % get element of the support
            a = master(mastertop(i,1),:);
            b = master(mastertop(i,3),:);
            % get shape function values and interpolation points coordinate
            slave_elems = find(elem_connectivity(i,:));
            % loop trough connected slave elements
            for j = slave_elems
                % get nodes
                c = slave(slavetop(j,1),:);
                d = slave(slavetop(j,3),:);
                % compute basis functions on master GP
                gSlave = ref2nod(gpRef,c,d);
                basisSlave = computeQuadraticSF(gpRef);
                % get Gauss Points in master reference coordinates
                gMasterRef = nod2ref(gSlave, a, b);
                basisMaster = computeQuadraticSF(gMasterRef);
                % check master GP that are not lying in the slave
                id = any([gMasterRef < -1-tol gMasterRef > 1+tol], 2);
                basisMaster(id,:) = 0;
                basisSlave = basisSlave.*gpWeight;
                matVal = basisMaster'*basisSlave;
                matVal = matVal*0.5*abs(c(1)-d(1));
                idMaster = [mastertop(i,1); mastertop(i,2); mastertop(i,3)];
                idSlave = [slavetop(j,1); slavetop(j,2); slavetop(j,3)];
                M_EB(idMaster, idSlave) = M_EB(idMaster, idSlave) + matVal;
            end
        end
    end


    %% COMPUTING MORTAR OPERATOR E = D^-1 * M
    E_RBF = D\M_RBF';
    E_EB = D\M_EB';
    E_EX = D\M_EX';

    %% INTERPOLATION TEST AND ERROR COMPUTATION
    % get lenght belonging to each slave node
    lNod = zeros(size(slave,1),1);
    for el = 1:size(slavetop,1)
        % get length of the element
        n1 = slavetop(el,1);
        n2 = slavetop(el,2);
        n3 = slavetop(el,3);
        l = sqrt((slave(n1,1) - slave(n3,1))^2+(slave(n1,2) - slave(n3,2))^2);
        lNod(n1) = lNod(n1) + l/4;
        lNod(n2) = lNod(n2) + l/2;
        lNod(n3) = lNod(n3) + l/4;
    end

    f = @(x) sin(4*x);
    fMaster = f(master(:,1)); % analytical function computed on master mesh
    plotSlave = (linspace(x1,x2,100))';
    fplotSlave = f(plotSlave);
    fSlave = f(slave(:,1)); % analytical function computed on master mesh

    errNormRBF(sizeCount) = computeInterpError(E_RBF,fMaster,fSlave,lNod);
    errNormEB(sizeCount) = computeInterpError(E_EB,fMaster,fSlave,lNod);
    errNormEX(sizeCount) = computeInterpError(E_EX,fMaster,fSlave,lNod);
end

%%

% limit_RBF = mean(sum(E_RBF_scale,2)-1);

% plot convergence of the interpolation
h = 1./(nMnodes*(2.^(0:nSizes-1))-1);
loglog(h,errNormEX,'r-o')
hold on
loglog(h,errNormEB,'g-o')
loglog(h,errNormRBF,'k-^')
legend('SB', 'EB', 'RBF')


% 
% % store results in mat-file
% if ~isfile("Results.mat")
%     outStruct= struct('nGP', nGP, 'nInt', nInt,...
%         'MeshSize', h, 'InterpolationError', errNormRBF);
%     save("Results.mat", "outStruct");
% else
%     out = load('Results.mat', "outStruct");
%     outStruct = out.outStruct; 
%     newStruct = struct('nGP', nGP, 'nInt', nInt,...
%         'MeshSize', h, 'InterpolationError', errNormRBF);
%     outStruct = [outStruct; newStruct];
%     save("Results.mat", "outStruct");
% end




%% PLOTTING INTERPOLATED FUNCTIONS using MAT FILE
% if all(master(:,2)==0)
%     c = 0.1;
% end
% figure(3)
% plot(master(:,1), c + master(:,2), 'b.-', 'LineWidth', 1.5,'MarkerSize',18) % master grid
% hold on 
% plot(slave(:,1), slave(:,2), 'r.-', 'LineWidth', 1.5, 'MarkerSize',18)
% ylim([-2 2])
% xlim([-2 2])
% plot(master, 3+fMaster, 'b')
% figure(2)
% plot(plotSlave, fplotSlave, 'k-', 'LineWidth', 1.2)
% hold on
% %plot(slave(:,1), fEX, 'k--')
% plot(slave(:,1), E_EX*fMaster, 'bs', 'MarkerSize', 10, 'MarkerFaceColor','b')
% if all(master(:,2)==0)
% plot(slave(:,1), E_EB*fMaster, 'r^', 'MarkerSize', 10)
% plot(slave(:,1), E_RBF*fMaster, 'gs', 'MarkerSize', 10)
% end
% legend('Exact function', 'SB integration', 'EB integration','RBF integration','Location','southeast')
% 
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 16);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% % export figure with quality
% % stmp = strcat('Images\', 'Intersection_interpolation', '.png');
% % exportgraphics(gcf,stmp,'Resolution',400)
% 
% figure(3)
% plot(nGP, errNormRBF, 'bo-')
% hold on
% plot(nGP, errNormEB, 'ro-')
% xlabel('Number of Gauss Points')
% ylabel('Norm of interpolation error')
% legend('RBF integration', 'Element-based integration')
% set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 16);
% a = get(gca,'XTickLabel');
% % set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% % % export figure with quality
% % stmp = strcat('Images\', 'err_vs_gp', '.png');
% % exportgraphics(gcf,stmp,'Resolution',400)
