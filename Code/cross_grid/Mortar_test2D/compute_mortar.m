function E = compute_mortar(masterInt, slaveInt, gauss, tagMaster, tagSlave)
% compute mortar matrix given from a pair of interfaces.
% INPUT: mesh objects, gauss class, edge tag for master and slave
% interfaces


% number of interpolation points for local master basis functions
nInt = 30;

% Gauss integration parameters
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
nNMaster = length(nodesMaster);
nodesSlave = unique(slaveInt.edges(slaveInt.edgeTag == tagSlave,:));
nNSlave = length(nodesSlave);
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

%% Compute slave matrix D 
for sID = 1:size(slavetop,1)
    % elem length
    i1 = slavetop(sID,1);
    i2 = slavetop(sID,2);
    h = norm(slave(i1,:)-slave(i2,:));
    Dloc = (h/6)*[2 1; 1 2];
    D([i1 i2], [i1 i2]) = D([i1 i2], [i1 i2]) + Dloc;
end
% extract only mass matrix entries belonging to the interface
D = D(nodesSlave, nodesSlave);

%% Approximate cross matrix computation

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
        basisF = 0.5 + gpRef*[-0.5 0.5];
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
M = M(nodesMaster, nodesSlave);

invD = inv(D);
E = invD*M';
end


