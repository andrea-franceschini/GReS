function [K,areaNod] = stiffPoisson(mesh, elem)
% computing classic stiffness matrix for Poisson problem
% uses element class from GReS to compute derivatives of basis functions
areaNod = zeros(mesh.nNodes,1);
nEntries = [9 16];
nSurfByType = histc(mesh.surfaceVTKType,[5, 9]);
iiVec = zeros(nEntries*nSurfByType,1);
jjVec = zeros(nEntries*nSurfByType,1);
KVec = zeros(nEntries*nSurfByType,1);
%
l1 = 0;
for el = 1:mesh.nSurfaces
    % Get the right material stiffness for each element
    switch mesh.surfaceVTKType(el)
        case 5 % Triangle
            N = getDerBasisF(elem.tri,el);
            vol = findArea(elem.tri,el);
            KLoc = N'*N*vol;
            s1 = 9;
    end
    % get global DoF
    % node indices
    nodes = mesh.surfaces(el,:);
    % update vector of nodal areas
    areaNod(nodes) = areaNod(nodes)+vol/3;
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

