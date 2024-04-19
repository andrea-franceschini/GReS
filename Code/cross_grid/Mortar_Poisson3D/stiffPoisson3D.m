function [K,volNod] = stiffPoisson3D(mesh, elem)
% computing classic stiffness matrix for Poisson problem
% uses element class from GReS to compute derivatives of basis functions
% for now the method works for hexahedrons only
volNod = zeros(mesh.nNodes,1);
nEntries = [16 64];
nCellByType = histc(mesh.cellVTKType,[10, 12]);
iiVec = zeros(nEntries*nCellByType,1);
jjVec = zeros(nEntries*nCellByType,1);
KVec = zeros(nEntries*nCellByType,1);
%
l1 = 0;
s1 = 64;
for el = 1:mesh.nCells
    % Get the right material stiffness for each element
    switch mesh.cellVTKType(el)
        case 12 % Hexahedron
            [N,dJWeighed] = elem.hexa.getDerBasisFAndDet(el,1);
            Ks = pagemtimes(N,'ctranspose',N,'none');
            Ks = Ks.*reshape(dJWeighed,1,1,[]);
            KLoc = sum(Ks,3);
            clear Ks;
    end
    % get global DoF
    % node indices
    nodes = mesh.cells(el,:);
    % update vector of nodal areas
    volNod(nodes) = volNod(nodes) + elem.hexa.findNodeVolume(el);
    % dof corresponding to dofs
    [jjLoc,iiLoc] = meshgrid(nodes',nodes');
    iiVec(l1+1:l1+s1) = iiLoc(:);
    jjVec(l1+1:l1+s1) = jjLoc(:);
    KVec(l1+1:l1+s1) = KLoc(:);
    l1 = l1 + s1;
end
% populate stiffness matrix in sparse format
K = sparse(iiVec, jjVec, KVec);
end

