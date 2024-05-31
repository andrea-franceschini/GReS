function K = computeStiffMat(mesh,element,D,GaussPts)
l1 = 0;
l2 = 0;
volNod = zeros(mesh.nNodes,1);
nEntryKLoc = (mesh.nDim^2)*(element.nNodesElem).^2;
for el=1:mesh.nCells
    % Get the right material stiffness for each element
    switch mesh.cellVTKType(el)
        case 10 % Tetrahedra
            N = getDerBasisF(element.tetra,el);
            vol = findVolume(element.tetra,el);
            B = zeros(6,4*mesh.nDim);
            B(element.indB(1:36,2)) = N(element.indB(1:36,1));
            % D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
            KLoc = B'*D*B*vol;
            %nodes = mesh.cells(el,:);
            %volNod(nodes) = volNod(nodes) + vol/4;
            s1 = nEntryKLoc(1);
            %fLoc = (B')*sz'*vol;
            s2 = 1;
            %             end
        case 12 % Hexahedra
            [N,dJWeighed] = getDerBasisFAndDet(element.hexa,el,1);
            B = zeros(6,8*mesh.nDim,GaussPts.nNode);
            B(element.indB(:,2)) = N(element.indB(:,1));
            Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
            Ks = Ks.*reshape(dJWeighed,1,1,[]);
            KLoc = sum(Ks,3);
            %nodes = mesh.cells(el,:);
            %volNod(nodes) = volNod(nodes) + element.hexa.findNodeVolume(el);
            clear Ks;
            s1 = nEntryKLoc(2);
            % fTmp = pagemtimes(B,'ctranspose',sz,'none');
            % fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
            % fLoc = sum(fTmp,3);
            s2 = GaussPts.nNode;
    end
    % get global DoF
    %dof = obj.dofm.ent2field('Poro',mesh.cells(el,1:mesh.cellNumVerts(el)));
    dof = getDoFID(mesh,el);
    [jjLoc,iiLoc] = meshgrid(dof,dof);
    iiVec(l1+1:l1+s1) = iiLoc(:);
    jjVec(l1+1:l1+s1) = jjLoc(:);
    KVec(l1+1:l1+s1) = KLoc(:);
    l1 = l1 + s1;
    l2 = l2 + s2;
end
% populate stiffness matrix
K = sparse(iiVec, jjVec, KVec);
K = full(K);
end


