function [E, M, D] = compute_mortar(masterInt, slaveInt, dirNod, nInt, nGP, tagMaster, tagSlave, tagMethod)
% compute mortar matrix given from a pair of interfaces.
% INPUT: mesh objects, gauss class, edge tag for master and slave
% interfaces

tolInt = 1.e-3;

% Gauss integration parameters
gauss = Gauss(12,nGP,1);
gpRef = gauss.coord;
gpWeight = gauss.weight;

% nodes coordinates
master = masterInt.coordinates(:,1:2);
slave = slaveInt.coordinates(:,1:2);

% edges topology
mastertop = masterInt.edges(masterInt.edgeTag == tagMaster, :);
slavetop = slaveInt.edges(slaveInt.edgeTag == tagSlave, :);


% number of slave and master nodes
nodesMaster = unique(masterInt.edges(masterInt.edgeTag == tagMaster,:));
nodesSlave = unique(slaveInt.edges(slaveInt.edgeTag == tagSlave,:));
% mass matrix on slave side
D = zeros(slaveInt.nNodes, slaveInt.nNodes);

% initializing mortar matrix resulting from RBF computation
M = zeros(masterInt.nNodes, slaveInt.nNodes);


% Element conntectivity matrix (contact search algorithm, that in the 1D interface case is basically trivial)
elem_connectivity = zeros(size(mastertop,1), size(slavetop,1));
for i = 1:size(mastertop,1)
    tmp = sort([master(mastertop(i,1),1),master(mastertop(i,2),1)]);
    a = tmp(1); b = tmp(2);
    % loop trough master element to find connectivity
    for j = 1:size(slavetop,1)
        tmp = sort([slave(slavetop(j,1),1),slave(slavetop(j,2),1)]);
        c = tmp(1); d = tmp(2);
        if ~any([a>d,c>b])
            % intersecting
            elem_connectivity(i,j) = 1;
        end
    end
end

% get edges with Dirichlet nodes (constant Basis Functions)
edgeDir = find(any(ismember(slavetop,dirNod),2));

%% Compute slave matrix D
for sID = 1:size(slavetop,1)
    % elem length
    i1 = slavetop(sID,1);
    i2 = slavetop(sID,2);
    h = norm(slave(i1,:)-slave(i2,:));
    if any(ismember(edgeDir,sID))
        Dloc = (h/2)*[1 1; 1 1];
    else
        Dloc = (h/6)*[2 1; 1 2];
    end
    D([i1 i2], [i1 i2]) = D([i1 i2], [i1 i2]) + Dloc;
end
% extract only mass matrix entries belonging to the interface
D = D(nodesSlave, nodesSlave);

%% Approximate cross matrix computation

switch tagMethod
    case "RBF_node" 
        % rbf interpolation is applied to the entire support of a master nodes
        for i = nodesMaster'
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
                if any(ismember(edgeDir,j))
                    basisF = ones(length(gpRef),2);
                else
                    basisF = 0.5 + gpRef.*[-0.5 0.5];
                end
                % product of the 2 functions in the element
                prod = basisF.*fSlave;
                % compute integrals with Gauss rule(use formula of basis change from ferro book)
                intVal = sum(prod.*gpWeight,1);
                intVal = intVal.*(0.5*sqrt((b(1)-a(1))^2 + (b(2)-a(2))^2));
                % populate cross matrix
                M(i,slavetop(j,1)) = M(i,slavetop(j,1)) + intVal(1);
                M(i,slavetop(j,2)) = M(i,slavetop(j,2)) + intVal(2);
            end
        end
    case "RBF"
        % rbf interpolation is applied to a single element for each of its
        % nodes
        for i = 1:size(mastertop,1)
            % get shape function values and interpolation points
            % coordinates
            [vals, ptsInt] = computeBasisF1D(i, nInt, mastertop, master);
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
                NMaster(isnan(NMaster)) = 0;
                NMaster(NMaster < 0) = 0;
                NMaster(NMaster > 1) = 0;
                % get basis function on Gauss points
                if any(ismember(edgeDir,j))
                    NSlave = ones(length(gpRef),2);
                else
                    NSlave = 0.5 + gpRef.*[-0.5 0.5];
                end
                % Compute local M matrix using Gauss integration
                Mloc = NMaster'*(NSlave.*gpWeight);
                % Apply the determinant of the jacobian (lenght of the
                % segment)
                Mloc = Mloc*(0.5*sqrt((b(1)-a(1))^2 + (b(2)-a(2))^2));
                % populate cross matrix
                idMaster = [mastertop(i,1); mastertop(i,2)];
                idSlave = [slavetop(j,1); slavetop(j,2)];
                M(idMaster, idSlave) = M(idMaster, idSlave) + Mloc;
            end
        end
    case "SB"
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
                if any(ismember(edgeDir,j))
                    basisSlave = ones(length(gSlave),2);
                else
                    basisSlave = 0.5 + gSlave.*[-0.5 0.5];
                end
                basisSlave = basisSlave.*gpWeight;
                matVal = basisMaster'*basisSlave;
                matVal = matVal.*(0.5*(i2(1)-i1(1)));
                idMaster = [mastertop(i,1); mastertop(i,2)];
                idSlave = [slavetop(j,1); slavetop(j,2)];
                M(idMaster, idSlave) = M(idMaster, idSlave) + matVal;
            end
        end
    case "EB"
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
                if any(ismember(edgeDir,j))
                    basisSlave = ones(length(gpRef),2);
                else
                    basisSlave = 0.5 + gpRef.*[-0.5 0.5];
                end
                % get Gauss Points in master reference coordinates
                gMasterRef = nod2ref(gSlave, a, b);
                basisMaster = 0.5 + gMasterRef*[-0.5 0.5];
                % check master GP that are not lying in the slave
                id = any([gMasterRef < -1-tolInt gMasterRef > 1+tolInt],2);
                basisMaster(id,:) = 0;
                basisSlave = basisSlave.*gpWeight;
                matVal = basisMaster'*basisSlave;
                matVal = matVal.*abs(0.5*(c(1)-d(1)));
                idMaster = [mastertop(i,1); mastertop(i,2)];
                idSlave = [slavetop(j,1); slavetop(j,2)];
                M(idMaster, idSlave) = M(idMaster, idSlave) + matVal;
            end
        end
end
%
M = M(nodesMaster, nodesSlave);
M = M';
E = D\M;
if strcmp(tagMethod,'RBF_node') 
    % rescaling the rows of the Mortar operator to improve convergence
    % in the case of node-wise interpolation
    E = E./sum(E,2);
end
end


