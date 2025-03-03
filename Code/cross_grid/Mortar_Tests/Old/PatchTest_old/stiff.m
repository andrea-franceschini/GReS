function K = stiff(mesh, elem, D, gauss)
% compute stiffness matrix for a certain domain
nEntryKLoc = [36 64];
nSurfByType = histc(mesh.surfaceVTKType,[5, 9]);
iiVec = zeros(nEntryKLoc*nSurfByType,1);
jjVec = zeros(nEntryKLoc*nSurfByType,1);
KVec = zeros(nEntryKLoc*nSurfByType,1);
%
l1 = 0;
l2 = 0;
for el = 1:mesh.nSurfaces
    % Get the right material stiffness for each element
    switch mesh.surfaceVTKType(el)
        case 5 % Triangle
            N = getDerBasisF(elem.tri,el);
            vol = findArea(elem.tri,el);
            B = zeros(3,6);
            B(elem.indB2D(1:12,2)) = N(elem.indB2D(1:12,1));
            KLoc = B'*D*B*vol;
            s1 = nEntryKLoc(1);
        case 9 % Hexahedra % CHECK LATER!
            [N,dJWeighed] = getDerBasisFAndDet(elem.quad,el,1);
            B = zeros(3,8,gauss.nNode);
            B(elem.indB2D(:,2)) = N(elem.indB2D(:,1));
            Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
            Ks = Ks.*reshape(dJWeighed,1,1,[]);
            KLoc = sum(Ks,3);
            s1 = nEntryKLoc(2);
    end
    % get global DoF
    % node indices
    nodes = mesh.surfaces(el,:);
    % dof corresponding to dofs
    dof = DofMap.getCompDoF(nodes);
    [jjLoc,iiLoc] = meshgrid(dof',dof');
    iiVec(l1+1:l1+s1) = iiLoc(:);
    jjVec(l1+1:l1+s1) = jjLoc(:);
    KVec(l1+1:l1+s1) = KLoc(:);
    l1 = l1 + s1;
    %l2 = l2 + s2;
end
% populate stiffness matrix in sparse format
K = sparse(iiVec, jjVec, KVec);
end


