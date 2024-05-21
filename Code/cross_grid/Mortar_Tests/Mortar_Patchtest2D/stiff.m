function K = stiff(mesh, elem, D)
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
            s2 = 1;
        case 12 % Hexahedra % CHECK LATER!
            % [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
            % B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
            % B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
            % [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
            %     state.conv.stress(l2+1:l2+obj.GaussPts.nNode,:), ...
            %     state.curr.strain(l2+1:l2+obj.GaussPts.nNode,:), ...
            %     dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:), el, state.t);
            % state.curr.status(l2+1:l2+obj.GaussPts.nNode,:) = status;
            % state.curr.stress((l2+1):(l2+obj.GaussPts.nNode),:) = sigma;
            % Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
            % Ks = Ks.*reshape(dJWeighed,1,1,[]);
            % KLoc = sum(Ks,3);
            % clear Ks;
            % s1 = obj.nEntryKLoc(2);
            % sz = sigma - state.iniStress(l2+1:l2+obj.GaussPts.nNode,:);
            % sz = reshape(sz',6,1,obj.GaussPts.nNode);
            % % fTmp = pagemtimes(B,'ctranspose',sz,'none');
            % % fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
            % % fLoc = sum(fTmp,3);
            % s2 = obj.GaussPts.nNode;
    end
    % get global DoF
    % node indices
    nodes = mesh.surfaces(el,:);
    % dof corresponding to dofs
    dof = getDoF(nodes);
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


