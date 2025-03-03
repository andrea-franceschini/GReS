function [E, M, D] = compute_mortar3D(mshMaster, mshSlave, conn, nGP, nInt)
% compute mortar operator using RBF
% INPUT: 3D interfaces as Mesh objects
% conn: connectivity matrix (returned by proper contact search algorithm)
% nGP: Number of Gauss Points for RBF integration
% nInt: number of interpolation points for each element of the master basis
% function support.


% Gauss class for approximated RBF integration
gauss = Gauss(12,nGP,2);

% nodes coordinates
master = mshMaster.coordinates;
slave = mshSlave.coordinates;

% surfaces topology
mastertop = mshMaster.surfaces;
slavetop = mshSlave.surfaces;


% number of slave and master nodes
nodesMaster = unique(mastertop);
nodesSlave = unique(slavetop);

% mass matrix on slave side
D = zeros(mshSlave.nNodes, mshMaster.nNodes);

% initializing mortar matrix resulting from RBF computation
M = zeros(mshMaster.nNodes, mshSlave.nNodes);

% Gauss integration for slave mass matrix
gaussMass = Gauss(12,3,2);
elemSlave = Elements(mshSlave, gaussMass);
%% Compute slave matrix D
l1 = 0;
s1 = 16;
for el = 1:size(slavetop,1)
    N1 = elemSlave.quad.getBasisFinGPoints();
    dJWeighed = elemSlave.quad.getDerBasisFAndDet(el,3);
    DLoc = N1'*diag(dJWeighed)*N1;
    dof = slavetop(el,:);
    [jjLoc,iiLoc] = meshgrid(dof,dof);
    iiVec(l1+1:l1+s1) = iiLoc(:);
    jjVec(l1+1:l1+s1) = jjLoc(:);
    DVec(l1+1:l1+s1) = DLoc(:);
    l1 = l1 + s1;
end
% extract only mass matrix entries belonging to the interface
D = sparse(iiVec, jjVec, DVec);

%% RBF approximate cross matrix computation (3D case)

elemMaster = Elements(mshMaster, gaussMass);

for i = nodesMaster'
    % get element of the support
    [master_elems,~] = find(mastertop == i);
    master_elems = sort(master_elems);
    % get shape function values and interpolation points coordinate
    [fMaster, ptsInt] = evalSF_3D(i, master_elems, nInt, mastertop, master, elemMaster);
    fiMM = zeros(length(ptsInt), length(ptsInt));
    % compute RBF weight of interpolant
    % local radius for radial basis function interpolation
    r = sqrt((max(ptsInt(:,1)) - min(ptsInt(:,1)))^2 + (max(ptsInt(:,2)) - min(ptsInt(:,2)))^2 + (max(ptsInt(:,3)) - min(ptsInt(:,3)))^2);
    for ii = 1:length(ptsInt)
        d = sqrt((ptsInt(:,1) - ptsInt(ii,1)).^2 + ((ptsInt(:,2)) - ptsInt(ii,2)).^2 + ((ptsInt(:,3)) - ptsInt(ii,3)).^2);
        rbf = pos(1-d./r).^4.*(1+d./r);
        fiMM(ii,:) = rbf;
    end
    % solve local system to get weight of interpolant
    weightF = fiMM\fMaster;
    weight1 = fiMM\ones(size(ptsInt,1),1);
    % get elements on slave surface connected to master support of node i
    slave_elems = [];
    for s = master_elems'
        slave_elems = [slave_elems; (find(conn(s,:)))'];
    end
    slave_elems = unique(slave_elems(:));

    % create istance of slave element class (for RBF integration)
    elemRBF = Elements(mshSlave, gauss);

    % loop trough connected slave elements
    for j = slave_elems'
        % get Gauss Points position in the real space
        ptsGauss = getGPointsLocation(elemRBF.quad,j);
        fiNM = zeros(size(ptsGauss,1), size(ptsInt,1));
        % compute interpolant on local Gauss points
        for jj = 1:size(ptsGauss,1)
            d = sqrt((ptsInt(:,1) - ptsGauss(jj,1)).^2 + ((ptsInt(:,2)) - ptsGauss(jj,2)).^2 + (ptsInt(:,3) - ptsGauss(jj,3)).^2);
            rbf = pos(1-d./r).^4.*(1+d./r);
            fiNM(jj,:) = rbf;
        end
        fMaster = (fiNM*weightF)./(fiNM*weight1);
        fMaster(isnan(fMaster)) = 0;
        fMaster(fMaster < 0) = 0;
        % get basis function on Gauss points
        basisF = getBasisFinGPoints(elemRBF.quad);
        % product of the 2 functions in the element
        prod = basisF.*fMaster;
        % compute integrals with Gauss rule(use formula of basis change from ferro book)
        dJWeighed = elemRBF.quad.getDerBasisFAndDet(j,3);
        intVal = sum(prod.*dJWeighed',1);
        % populate cross matrix
        M(i,slavetop(j,:)) = M(i,slavetop(j,:)) + intVal;
    end
end
M = M';
%
E = D\M;
end


