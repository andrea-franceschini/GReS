function E = compute_mortar_SB(masterInt, slaveInt, tagMaster, tagSlave, nGP)
% compute mortar matrix given from a pair of interfaces.
% INPUT: mesh objects, gauss class, edge tag for master and slave
% interfaces

% uses Segment based approach for exact computation of the mortar integral
% valid only if the connected interfaces are flat

tol = 1e-3;

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


% Exact cross matrix computation (DO ONLY IN THE 1D CASE)

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
        matVal = matVal.*(0.5*(i2(1)-i1(1)));
        idMaster = [mastertop(i,1); mastertop(i,2)];
        idSlave = [slavetop(j,1); slavetop(j,2)];
        M(idMaster, idSlave) = M(idMaster, idSlave) + matVal;
    end
end

M = M(nodesMaster, nodesSlave);

invD = inv(D);
E = invD*M';
end


