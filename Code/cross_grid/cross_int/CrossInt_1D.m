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

%% INPUT DATAS

close all; clear

% INPUT
degree = 1;
tol = 1.e-4;

% master and slave surfaces node position along X axis

nSizes = 1;
errNormRBF_scale = zeros(nSizes,1);
nMnodes = 12; % number of nodes for the first mesh refinment
errNormEX = errNormRBF_scale;
errNormEB = errNormRBF_scale;
errNormRBF_noScale = errNormRBF_scale;
errNormRBF_elem_noScale = errNormRBF_scale;
errNormRBF_elem = errNormRBF_scale;
for sizeCount = 1:nSizes
    nMaster = nMnodes*(2^(sizeCount-1));
    nSlave = round((2/3)*nMaster);
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
    fact = 0.5;
    master(:,2) = -fact*sm*rand(nMaster,1)+fact*sm*rand(nMaster,1);
    slave(:,2) = -fact*ss*rand(nSlave,1)+fact*ss*rand(nSlave,1);


    % number of RBF interpolation points for each element
    nInt = 5;

    % Number of integration points for RBF testing (GP class taken from GReS)
    nGP = 4;

    % Inizialize output matrices
    D = zeros(length(slave), length(slave));

    % initializing mortar matrix resulting from RBF computation

    M_RBF = zeros(length(master), length(slave));
    M_RBF_elem = zeros(length(master), length(slave));
    M_EB = zeros(length(master), length(slave));
    M_EX = zeros(size(master,1), size(slave,1));
    % Build a topology matrix for master/slave surfs based on nodes position
    mastertop = build_topol(master(:,1));
    slavetop = build_topol(slave(:,1));

    % Element conntectivity matrix
    elem_connectivity = zeros(size(mastertop,1), size(slavetop,1));
    for i = 1:size(mastertop,1)
        a = master(mastertop(i,1));
        b = master(mastertop(i,2));
        % loop trough master element to find connectivity
        for j = 1:size(slavetop,1)
            c = slave(slavetop(j,1));
            d = slave(slavetop(j,2));
            if ~any([a>d,c>b])
                % intersecting
                elem_connectivity(i,j) = 1;
            end
        end
    end

    %% COMPUTE SLAVE MATRIX D
    for sID = 1:size(slavetop,1)
        % elem length
        i1 = slavetop(sID,1);
        i2 = slavetop(sID,2);
        h = norm(slave(i1,:)-slave(i2,:));
        Dloc = (h/6)*[2 1; 1 2];
        D([i1 i2], [i1 i2]) = D([i1 i2], [i1 i2]) + Dloc;
    end

    %% Node wise RBF CROSS MATRIX COMPUTATION
    tic
    for i = 1:nMaster
        % get element of the support
        [master_elems,~] = find(mastertop == i);
        master_elems = sort(master_elems);
        % get shape function values and interpolation points coordinate
        [fMaster, ptsInt] = evalSF(i, master_elems, nInt, mastertop, master);
        fiMM = zeros(length(ptsInt), length(ptsInt));
        % compute RBF weight of interpolant
        % local radius for radial basis function interpolation
        r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2);
        for ii = 1:length(ptsInt)
            d = sqrt((ptsInt(:,1) - ptsInt(ii,1)).^2 + ((ptsInt(:,2)) - ptsInt(ii,2)).^2);
            rbf = pos(1-d./r).^4.*(1+d./r);
            fiMM(ii,:) = rbf;
        end
        % solve local system to get weight of interpolant
        weightF = fiMM\fMaster;
        weight1 = fiMM\ones(size(ptsInt,1),1);
        % get elements on slave surface connected to master support of node i
        slave_elems = [];
        for s = master_elems'
            slave_elems = [slave_elems; (find(elem_connectivity(s,:)))'];
        end
        slave_elems = unique(slave_elems(:));
        % loop trough connected slave elements
        for j = slave_elems'
            % get nodes of slave elements
            a = slave(slavetop(j,1),:);
            b = slave(slavetop(j,2),:);
            for gc = 1:length(nGP)
                g = Gauss(12,nGP,1);
                gpRef = g.coord;
                gpWeight = g.weight;
                % get Gauss Points position in the real space
                ptsGauss = ref2nod(gpRef, a, b);
                fiNM = zeros(size(ptsGauss,1), size(ptsInt,1));
                % compute interpolant on local Gauss points
                for jj = 1:size(ptsGauss,1)
                    d = sqrt((ptsInt(:,1) - ptsGauss(jj,1)).^2 + ((ptsInt(:,2)) - ptsGauss(jj,2)).^2);
                    rbf = pos(1-d./r).^4.*(1+d./r);
                    fiNM(jj,:) = rbf;
                end
                fSlave = (fiNM*weightF)./(fiNM*weight1);
                fSlave(isnan(fSlave)) = 0;
                fSlave(fSlave < 0) = 0;
                % get basis function on Gauss points
                basisF = 0.5 + gpRef*[-0.5 0.5];
                % product of the 2 functions in the element
                prod = basisF.*fSlave;
                % compute integrals with Gauss rule(use formula of basis change from ferro book)
                intVal = sum(prod.*gpWeight,1);
                intVal = intVal.*(0.5*sqrt((b(1)-a(1))^2 + (b(2)-a(2))^2));
                % populate cross matrix
                M_RBF(i,slavetop(j,1)) = M_RBF(i,slavetop(j,1)) + intVal(1);
                M_RBF(i,slavetop(j,2)) = M_RBF(i,slavetop(j,2)) + intVal(2);
            end
        end
    end

    %% Element-wise RBF matrix computation  
    % rbf interpolation is applied to a single element for each of its
    % nodes
    % It's more efficient for different reasons:
    % a smoother function is interpolated (no kinks in the node)
    % smaller linear systems are computed
    % for the rescaling, only one system is solved for each element
    for i = 1:size(mastertop,1)
        % get shape function values and interpolation points
        % coordinates
        [vals, ptsInt] = computeBasisF1D(i, 2*nInt, mastertop, master, degree);
        % compute interpolation matrix
        fiMM = zeros(length(ptsInt), length(ptsInt));
        % compute RBF weight of interpolant
        % circumradius of the master element
        r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2);
        for ii = 1:length(ptsInt)
            d = sqrt((ptsInt(:,1) - ptsInt(ii,1)).^2 + ((ptsInt(:,2)) - ptsInt(ii,2)).^2);
            rbf = pos(1-d./r).^4.*(1+d./r);
            fiMM(ii,:) = rbf;
        end
        % solve local system to get weight of interpolant
        weightF = fiMM\vals;
        weight1 = fiMM\ones(size(ptsInt,1),1);
        % get slave elements connected to the master
        slave_elems = find(elem_connectivity(i,:));
        % loop trough connected slave elements
        for j = slave_elems
            % get nodes of slave elements
            a = slave(slavetop(j,1),:);
            b = slave(slavetop(j,2),:);
            % get Gauss Points position in the real space
            ptsGauss = ref2nod(gpRef, a, b);
            fiNM = zeros(size(ptsGauss,1), size(ptsInt,1));
            % compute interpolant on local Gauss points
            for jj = 1:size(ptsGauss,1)
                d = sqrt((ptsInt(:,1) - ptsGauss(jj,1)).^2 + ((ptsInt(:,2)) - ptsGauss(jj,2)).^2);
                rbf = pos(1-d./r).^4.*(1+d./r);
                fiNM(jj,:) = rbf;
            end
            NMaster = (fiNM*weightF)./(fiNM*weight1);
            %NMaster = fiNM*weightF;
            % automatically detect GP outside the master support
            NMaster(isnan(NMaster)) = 0;
            NMaster(NMaster < 0) = 0;
            NMaster(NMaster > 1) = 0;
            % get basis function on Gauss points
            NSlave = 0.5 + gpRef.*[-0.5 0.5];
            % Compute local M matrix using Gauss integration
            Mloc = NMaster'*(NSlave.*gpWeight);
            % Apply the determinant of the jacobian (lenght of the
            % segment)
            Mloc = Mloc*(0.5*sqrt((b(1)-a(1))^2 + (b(2)-a(2))^2));
            % populate cross matrix
            idMaster = [mastertop(i,1); mastertop(i,2)];
            idSlave = [slavetop(j,1); slavetop(j,2)];
            M_RBF_elem(idMaster, idSlave) = M_RBF_elem(idMaster, idSlave) + Mloc;
        end
    end
    nodT = toc;

    %% SEGMENT-BASED CROSS MATRIX COMPUTATION (valid for flat interfaces)
    tic
    if all(master(:,2)==0)
        gEx = Gauss(12,3,1);
        gpRef = gEx.coord;
        gpWeight = gEx.weight;
        for i = 1:size(mastertop,1)
            % get element of the support
            a = master(mastertop(i,1),:);
            b = master(mastertop(i,2),:);
            % get shape function values and interpolation points coordinate
            slave_elems = find(elem_connectivity(i,:));
            % loop trough connected slave elements
            for j = slave_elems
                % get nodes
                c = slave(slavetop(j,1),:);
                d = slave(slavetop(j,2),:);
                % find intersection between master and slave (just looking at the
                % x-projection)
                tmp = sort([a(1) b(1) c(1) d(1)]);
                i1 = [tmp(2) 0];
                i2 = [tmp(3) 0];
                % get GP in the intersection element
                gPts = ref2nod(gpRef, i1, i2);
                % compute basis function on master element points
                gMaster = nod2ref(gPts, a, b);
                basisMaster = 0.5 + gMaster*[-0.5 0.5];
                % compute basis function on slave elements points
                gSlave = nod2ref(gPts, c, d);
                basisSlave = 0.5 + gSlave.*[-0.5 0.5];
                basisSlave = basisSlave.*gpWeight;
                matVal = basisMaster'*basisSlave;
                matVal = matVal*0.5*abs(i2(1)-i1(1));
                idMaster = [mastertop(i,1); mastertop(i,2)];
                idSlave = [slavetop(j,1); slavetop(j,2)];
                M_EX(idMaster, idSlave) = M_EX(idMaster, idSlave) + matVal;
            end
        end
    end
    elemT = toc;

    %% ELEMENT-BASED CROSS MATRIX COMPUTATION (valid for flat interfaces)

    if all(master(:,2)==0)
        for gc = 1:length(nGP)
            gEx = Gauss(12,nGP,1);
            gpRef = gEx.coord;
            gpWeight = gEx.weight;
            for i = 1:size(mastertop,1)
                % get element of the support
                a = master(mastertop(i,1),:);
                b = master(mastertop(i,2),:);
                % get shape function values and interpolation points coordinate
                slave_elems = find(elem_connectivity(i,:));
                % loop trough connected slave elements
                for j = slave_elems
                    % get nodes
                    c = slave(slavetop(j,1),:);
                    d = slave(slavetop(j,2),:);
                    % compute basis functions on master GP
                    gSlave = ref2nod(gpRef,c,d);
                    basisSlave = 0.5 + gpRef*[-0.5 0.5];
                    % get Gauss Points in master reference coordinates
                    gMasterRef = nod2ref(gSlave, a, b);
                    basisMaster = 0.5 + gMasterRef*[-0.5 0.5];
                    % check master GP that are not lying in the slave
                    id = any([gMasterRef < -1-tol gMasterRef > 1+tol], 2);
                    basisMaster(id,:) = 0;
                    basisSlave = basisSlave.*gpWeight;
                    matVal = basisMaster'*basisSlave;
                    matVal = matVal*0.5*abs(c(1)-d(1));
                    idMaster = [mastertop(i,1); mastertop(i,2)];
                    idSlave = [slavetop(j,1); slavetop(j,2)];
                    M_EB(idMaster, idSlave) = M_EB(idMaster, idSlave) + matVal;
                end
            end
        end
    end


    %% COMPUTING MORTAR OPERATOR E = D^-1 * M
    E_RBF_noScale = D\M_RBF';
    E_RBF_scale = E_RBF_noScale./sum(E_RBF_noScale,2);
    E_RBF_elem = D\M_RBF_elem';
    E_RBF_elem = E_RBF_elem./sum(E_RBF_elem,2);
    E_EB = D\M_EB';
    E_EX = D\M_EX';

    %% INTERPOLATION TEST AND ERROR COMPUTATION
    % get lenght belonging to each slave node
    lNod = zeros(size(slave,1),1);
    for el = 1:size(slavetop,1)
        % get length of the element
        n1 = slavetop(el,1);
        n2 = slavetop(el,2);
        l = sqrt((slave(n1,1) - slave(n2,1))^2+(slave(n1,2) - slave(n2,2))^2);
        lNod(n1) = lNod(n1) + l/2;
        lNod(n2) = lNod(n2) + l/2;
    end

    f = @(x) sin(3*x) + exp(2*x);
    fMaster = f(master(:,1)); % analytical function computed on master mesh
    plotSlave = (linspace(x1,x2,100))';
    fplotSlave = f(plotSlave);
    fSlave = f(slave(:,1)); % analytical function computed on master mesh

    errNormRBF_scale(sizeCount) = computeInterpError(E_RBF_scale,fMaster,fSlave,lNod);
    errNormRBF_noScale(sizeCount) = computeInterpError(E_RBF_noScale,fMaster,fSlave,lNod);
    [errNormRBF_elem(sizeCount), fRBF] = computeInterpError(E_RBF_elem,fMaster,fSlave,lNod);
    errNormEB(sizeCount) = computeInterpError(E_EB,fMaster,fSlave,lNod);
    errNormEX(sizeCount) = computeInterpError(E_EX,fMaster,fSlave,lNod);
end

fName = strcat('errRBFint',num2str(nInt));
f = fopen(fName,"w");
fprintf(f, '%2.5e \n',errNormRBF_elem);

%%

limit_RBF = mean(sum(E_RBF_scale,2)-1);
errRBF2 = load('errRBFint2');
errRBF4 = load('errRBFint4');
errRBF8 = load('errRBFint8');
% plot convergence of the interpolation
h = 1./(nMnodes*2.^(0:nSizes-1));
loglog(h,errNormEX,'r-o')
hold on
loglog(h,errNormEB,'g-o')
%loglog(h,errNormRBF_scale,'b-o')
%loglog(h,errNormRBF_noScale,'b--o')
loglog(h,errNormRBF_elem,'k-^')
legend('SB', 'EB', 'RBF - Elem')

% errRBF2 = load('errRBFint2');
% errRBF4 = load('errRBFint4');
% errRBF8 = load('errRBFint8');
% % plot convergence of the interpolation
% h = 1./(nMnodes*2.^(0:nSizes-1));
% %oglog(h,errNormEX,'r-o')
% %loglog(h,errNormEB,'g-o')
% loglog(h,errRBF2,'r-o')
% hold on
% loglog(h,errRBF4,'b-o')
% loglog(h,errRBF8,'k-o')
% %loglog(h,errNormRBF_noScale,'b--o')
% %loglog(h,errNormRBF_elem,'k-^')
% legend('nInt = 2','nInt = 4', 'nInt = 8')


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
c = 0;
if all(master(:,2)==0)
    c = 0.1;
end
figure(3)
plot(master(:,1), c + master(:,2), 'b.-', 'LineWidth', 1.5,'MarkerSize',18) % master grid
hold on 
plot(slave(:,1), slave(:,2), 'r.-', 'LineWidth', 1.5, 'MarkerSize',18)
ylim([-2 2])
xlim([-2 2])
%plot(master, fMaster, 'b')
figure(2)
plot(plotSlave, fplotSlave, 'k-', 'LineWidth', 1.2)
hold on
%plot(slave(:,1), fEX, 'k--')
plot(slave(:,1), fRBF, 'bs', 'MarkerSize', 10, 'MarkerFaceColor','b')
if all(master(:,2)==0)
plot(slave(:,1), fEX, 'r^', 'MarkerSize', 10)
plot(slave(:,1), fEB, 'gs', 'MarkerSize', 10)
end
legend('Exact function', 'RBF integration', 'Segment-based integration','Element-based integration','Location','southeast')

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
