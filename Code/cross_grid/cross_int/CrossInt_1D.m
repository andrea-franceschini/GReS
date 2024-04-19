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

tol = 1.e-2;

% master and slave surfaces node position along X axis
nMaster = 8;
nSlave = 8;
master = zeros(nMaster,2); % coordinates of master side
slave = zeros(nSlave,2);
master(:,1) = (linspace(-1,1,nMaster))';
slave(:,1) = (linspace(-1,1,nSlave))';

% y-axis (to add possible overlapping bewteen the interfaces)
% we set the maximum overlapping to 1/10 of the grid size
sm = abs(master(2)-master(1));
ss = abs(slave(2)-slave(1));
fact = 0;
master(:,2) = -fact*sm*rand(nMaster,1)+fact*sm*rand(nMaster,1);
slave(:,2) = -fact*ss*rand(nSlave,1)+fact*ss*rand(nSlave,1);


 % number of RBF interpolation points for each element
nInt = 100;

% Number of integration points for RBF testing (GP class taken from GReS)
nGP = 6;

% Inizialize output matrices
D = zeros(length(slave), length(slave));
M_RBF = cell(length(nGP),1);
E_RBF = cell(length(nGP),1);
fRBF = cell(length(nGP),1);
M_EB = cell(length(nGP),1);
E_EB = cell(length(nGP),1);
fEB = cell(length(nGP),1);
% initializing mortar matrix resulting from RBF computation
for c = 1:length(nGP)
    M_RBF{c} = zeros(length(master), length(slave));
    E_RBF{c} = zeros(length(slave), length(master));
    M_EB{c} = zeros(length(master), length(slave));
    E_EB{c} = zeros(length(slave), length(master));
end
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

%% RBF CROSS MATRIX COMPUTATION
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
            g = Gauss(12,nGP(gc),1);
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
            M_RBF{gc}(i,slavetop(j,1)) = M_RBF{gc}(i,slavetop(j,1)) + intVal(1);
            M_RBF{gc}(i,slavetop(j,2)) = M_RBF{gc}(i,slavetop(j,2)) + intVal(2);
        end
    end
end

%% SEGMENT-BASED CROSS MATRIX COMPUTATION (valid for flat interfaces)
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
            matVal = matVal.*abs(0.5*(i2(1)-i1(1)));
            idMaster = [mastertop(i,1); mastertop(i,2)];
            idSlave = [slavetop(j,1); slavetop(j,2)];
            M_EX(idMaster, idSlave) = M_EX(idMaster, idSlave) + matVal;
        end
    end
end

%% ELEMENT-BASED CROSS MATRIX COMPUTATION (valid for flat interfaces)

if all(master(:,2)==0)
    for gc = 1:length(nGP)
        gEx = Gauss(12,nGP(gc),1);
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
                matVal = matVal.*abs(0.5*(c(1)-d(1)));
                idMaster = [mastertop(i,1); mastertop(i,2)];
                idSlave = [slavetop(j,1); slavetop(j,2)];
                M_EB{gc}(idMaster, idSlave) = M_EB{gc}(idMaster, idSlave) + matVal;
            end
        end
    end
end


%% COMPUTING MORTAR OPERATOR E = D^-1 * M
invD = inv(D);
for c = 1:length(nGP)
    E_RBF{c} = invD*(M_RBF{c})';
    E_EB{c} = invD*(M_EB{c})';
end
E_EX = invD*(M_EX)';


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

f = @(x) 3 + 0.000001*x;
fMaster = f(master(:,1)); % analytical function computed on master mesh
plotSlave = (linspace(-1,1,100))';
fplotSlave = f(plotSlave);
fSlave = f(slave(:,1)); % analytical function computed on master mesh
errNormRBF = zeros(length(nGP),1);
errNormEB = zeros(length(nGP),1);
for c = 1:length(nGP)
    fRBF{c} = E_RBF{c} * fMaster;
    err2RBF = (fSlave - fRBF{c}).^2;
    errNormRBF(c) = sqrt(sum(err2RBF.*lNod));
    fEB{c} = E_EB{c} * fMaster;
    err2EB = (fSlave - fEB{c}).^2;
    errNormEB(c) = sqrt(sum(err2EB.*lNod));
end
fEX = E_EX * fMaster;
err2EX = (fSlave - fEX).^2;
errNormEX = sqrt(sum(err2EX.*lNod));


%% PLOTTING INTERPOLATED FUNCTIONS
figure(1)
plot(master(:,1), master(:,2), 'b.-', 'LineWidth', 1.5,'MarkerSize',18) % master grid
hold on 
plot(slave(:,1), 0.3 + slave(:,2), 'r.-', 'LineWidth', 1.5, 'MarkerSize',18)
ylim([-2 2])
xlim([-2 2])
% plot(master, 3+fMaster, 'b')
figure(2)
plot(plotSlave, fplotSlave, 'k-', 'LineWidth', 1.2)
hold on
%plot(slave(:,1), fEX, 'k--')
plot(slave(:,1), fRBF{end}, 'bs', 'MarkerSize', 10, 'MarkerFaceColor','b')
if all(master(:,2)==0)
plot(slave(:,1), fEX, 'r^', 'MarkerSize', 10)
plot(slave(:,1), fEB{end}, 'gs', 'MarkerSize', 10)
end
legend('Exact function', 'RBF integration', 'Segment-based integration','Element-based integration','Location','southeast')

set(findall(gcf, 'type', 'text'), 'FontName', 'Liberation Serif','FontSize', 16);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName', 'Liberation Serif','FontSize', 10)
% export figure with quality
stmp = strcat('Images\', 'Intersection_interpolation', '.png');
exportgraphics(gcf,stmp,'Resolution',400)

figure(3)
plot(nGP, errNormRBF, 'bo-')
hold on
plot(nGP, errNormEB, 'ro-')
xlabel('Number of Gauss Points')
ylabel('Norm of interpolation error')
legend('RBF integration', 'Element-based integration')









