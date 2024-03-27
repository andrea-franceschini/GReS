function E = compute_mortar_EB(masterInt, slaveInt, gauss, tagMaster, tagSlave)
% compute mortar matrix given from a pair of interfaces.
% INPUT: mesh objects, gauss class, edge tag for master and slave
% interfaces

% uses Segment based approach for exact computation of the mortar integral
% valid only if the connected interfaces are flat

tol = 1e-3;
tolInt = 1.e-3;

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
nodesSlave = unique(slaveInt.edges(slaveInt.edgeTag == tagSlave,:));
%
if abs(max(masterInt.coordinates(nodesMaster,2))-min(masterInt.coordinates(nodesMaster,2))) > tol
    return
end

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


% Element based cross matrix computation (DO ONLY IN THE 1D CASE)

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


M = M(nodesMaster, nodesSlave);

invD = inv(D);
E = invD*M';
end


