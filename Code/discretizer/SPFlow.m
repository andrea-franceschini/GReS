classdef SPFlow < handle
    %POROMECHANICS
    % Subclass of Discretizer
    % Implement Poromechanics methods to assemble the stiffness matrix and
    % the residual contribution

    properties
        model
        simParams
        dofm
        mesh
        elements
        faces
        material
        GaussPts
        trans
        isIntFaces
        H               % stiffness matrix
        P               % capacity matrix
        rhs             % residual
        rhsGrav         % gravity contribute to rhs
    end

    methods (Access = public)
        function obj = SPFlow(symmod,params,dofManager,grid,mat,data)
            obj.model = symmod;
            obj.simParams = params;
            obj.dofm = dofManager;
            obj.mesh = grid.topology;
            obj.elements = grid.cells;
            obj.faces = grid.faces;
            obj.material = mat;
            if ~isempty(data)
                obj.GaussPts = data{1};
            end

            if obj.model.isFVTPFABased('Flow')
                obj.computeTrans;
                %get cells with active flow model
                % internal faces inside each subdomain
                nSub = length(obj.dofm.subDomains);
                IntFaces = zeros(obj.faces.nFaces,nSub);
                flowCells = [];
                for i = 1:nSub
                    if any(strcmp(obj.dofm.subDomains(i).physics,"SPFlow"))
                        flowCells = [flowCells; find(obj.dofm.subCells(:,i))];
                    end
                end

                for i = 1:nSub
                    intcells = find(obj.dofm.subCells(:,i));
                    tmp1 = any(ismember(obj.faces.faceNeighbors, intcells), 2);
                    tmp2 = all(ismember(obj.faces.faceNeighbors, flowCells), 2);
                    IntFaces(:,i) = all([tmp1 tmp2],2);
                end
                obj.isIntFaces = logical(IntFaces);
            end

            computeRHSGravTerm(obj);

        end


        function computeMat(obj,varargin)
            if isempty(obj.H) && isempty(obj.P)
                if obj.model.isFEMBased('Flow')
                    computeMatFEM(obj,varargin{:});
                elseif obj.model.isFVTPFABased('Flow')
                    mu = obj.material.getFluid().getDynViscosity();
                    computeStiffMatFV(obj,1/mu);
                    computeCapMatFV(obj);
                end
            end
        end

        function computeMatFEM(obj,varargin) %provisional method exploiting dof manager workflow
            % dealing with input params
            if ~isempty(varargin)
               K = varargin{1};
               poro = varargin{2};
               alpha = varargin{3};
               if numel(poro)==1
                  porosity = repmat(poro,obj.mesh.nCells,1);
               end
               if numel(poro)==1
                  rockComp = repmat(alpha,obj.mesh.nCells,1);
               end
            end
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'SPFlow'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            %nSubCells = length(subCells); %number of cells in subdomain
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            % Compute the stiffness (H) and mass (P) matrices for the flow problem by FEM
            iiVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            jjVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            HVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            PVec = zeros((obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            % Get the fluid compressibility
            beta = obj.material.getFluid().getFluidCompressibility();
            if nSubCellsByType(2) > 0
                N1 = obj.elements.hexa.getBasisFinGPoints();
            end
            % Get the fluid dynamic viscosity
            mu = obj.material.getFluid().getDynViscosity();
            %
            l1 = 0;
            for el = subCells'
                % Get the rock permeability, porosity and compressibility 
                if isempty(varargin) % input from PorousRock class 
                   permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
                   poro = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPorosity();
                   if  any(strcmp(obj.dofm.subDomains(subInd).physics,'Poro'))
                      alpha = 0; %this term is not needed in coupled formulation
                   else
                      alpha = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
                      %solid skeleton contribution to storage term as oedometric compressibility .
                   end
                else % direct input
                   permMat = diag(repmat(K(el),3,1));
                   poro = porosity(el);
                   alpha = rockComp(el);
                end
                % Compute the element matrices based on the element type
                % (tetrahedra vs. hexahedra)
                switch obj.mesh.cellVTKType(el)
                    case 10 % Tetrahedra
                        % Computing the H matrix contribution
                        N = obj.elements.tetra.getDerBasisF(el);
                        %               vol = getVolume(obj.elements,el);
                        HLoc = N'*permMat*N*obj.elements.vol(el)/mu;
                        s1 = obj.elements.nNodesElem(1)^2;
                        % Computing the P matrix contribution
                        PLoc = ((alpha + poro*beta)*obj.elements.vol(el)/20)*(ones(obj.elements.nNodesElem(1))...
                            + eye(obj.elements.nNodesElem(1)));
                    case 12 % Hexa
                        [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                        permMat = permMat/mu;
                        Hs = pagemtimes(pagemtimes(N,'ctranspose',permMat,'none'),N);
                        Hs = Hs.*reshape(dJWeighed,1,1,[]);
                        HLoc = sum(Hs,3);
                        clear Hs;
                        s1 = obj.elements.nNodesElem(2)^2;
                        % Computing the P matrix contribution
                        PLoc = (alpha+poro*beta)*(N1'*diag(dJWeighed)*N1);
                end
                %Getting dof associated to Flow subphysic
                dof = obj.dofm.ent2field('SPFlow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                HVec(l1+1:l1+s1) = HLoc(:);
                PVec(l1+1:l1+s1) = PLoc(:);
                l1 = l1 + s1;
            end
            % Assemble H and P matrices defined as new fields of
            obj.H = sparse(iiVec, jjVec, HVec);
            obj.P = sparse(iiVec, jjVec, PVec);
        end


        function computeStiffMatFV(obj,lw)
            % Inspired by MRST
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'SPFlow'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCells = length(subCells); 
            %get pairs of faces that contribute to the subdomain
            intFaces = any(obj.isIntFaces(:,subInd),2);
            neigh1 = obj.faces.faceNeighbors(intFaces,1);
            neigh2 = obj.faces.faceNeighbors(intFaces,2);
            neigh1dof = obj.dofm.ent2field('SPFlow',neigh1);
            neigh2dof = obj.dofm.ent2field('SPFlow',neigh2);
            % Transmissibility of internal faces
            tmpVec = lw.*obj.trans(intFaces);
            % tmpVec = lw.*tmpVec;
            sumDiagTrans = accumarray([neigh1dof; neigh2dof], ...
                repmat(tmpVec,[2,1]),[nSubCells,1]);
            % Assemble H matrix
            obj.H = sparse([neigh1dof; neigh2dof; (1:nSubCells)'], ...
                [neigh2dof; neigh1dof; (1:nSubCells)'], ...
                [-tmpVec; -tmpVec; ...
                sumDiagTrans]);
        end



        function computeCapMatFV(obj,varargin)
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'SPFlow'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCells = length(subCells);
            poroMat = zeros(nSubCells,1);
            alphaMat = zeros(nSubCells,1);
            beta = obj.material.getFluid().getFluidCompressibility();
            for m = 1:obj.mesh.nCellTag
                if ~any(strcmp(obj.dofm.subDomains(subInd).physics,'Poro'))
                    % compute alpha only if there's no coupling in the
                    % subdomain
                    alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
                end
                poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
            end
            % (alpha+poro*beta)
            PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));
            PVal = PVal.*obj.elements.vol(subCells);
            dof = obj.dofm.ent2field('SPFlow', subCells);
            obj.P = sparse(dof,dof,PVal);
        end


        function computeRhs(obj,stateTmp,statek,dt)
            % Compute the residual of the flow problem
            lw = obj.material.getFluid().getDynViscosity();
            theta = obj.simParams.theta;
            ents = obj.dofm.field2ent('SPFlow');
            rhsStiff = theta*obj.H*stateTmp.pressure(ents) + (1-theta)*obj.H*statek.pressure(ents);
            rhsCap = (obj.P/dt)*(stateTmp.pressure(ents) - statek.pressure(ents));
            obj.rhs = rhsStiff + rhsCap;
            gamma = obj.material.getFluid().getFluidSpecWeight();
            %adding gravity rhs contribute
            if gamma > 0
                if isFEMBased(obj.model,'Flow')
                    obj.rhs = obj.rhs + obj.rhsGrav;
                elseif isFVTPFABased(obj.model,'Flow')
                    obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,lw);
                end
            end
        end

        function computeRHSGravTerm(obj)
            % Compute the gravity contribution
            % Get the fluid specific weight and viscosity'
            rhsTmp = zeros(sum(obj.dofm.numDof),1);
            gamma = obj.material.getFluid().getFluidSpecWeight();
            if gamma > 0
                subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'SPFlow'));
                [subCells, ~] = find(obj.dofm.subCells(:,subInd));
                if isFEMBased(obj.model,'Flow')
                    for el = subCells'
                        % Get the material permeability
                        permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
                        %             permMat = permMat/mu;
                        switch obj.mesh.cellVTKType(el)
                            case 10 % Tetrahedra
                                N = obj.elements.tetra.getDerBasisF(el);
                                fLoc = (N'*permMat(:,3))*obj.elements.vol(el)*gamma;
                            case 12 % Hexa
                                [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                                fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
                                fs = fs.*reshape(dJWeighed,1,1,[]);
                                fLoc = sum(fs,3)*gamma;
                        end
                        %
                        dof = obj.dofm.ent2field('SPFlow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                        rhsTmp(dof) = rhsTmp(dof) + fLoc; 
                    end
                    obj.rhsGrav = rhsTmp(obj.dofm.getDoF('SPFlow'));
                elseif isFVTPFABased(obj.model,'Flow')
                    intFaces = any(obj.isIntFaces(:,subInd),2);
                    neigh = obj.faces.faceNeighbors(intFaces,:);
                    zVec = obj.elements.cellCentroid(:,3);
                    zNeigh = zVec(neigh);
                    obj.rhsGrav = gamma*obj.trans(obj.isIntFaces(:,subInd)).*(zNeigh(:,1) - zNeigh(:,2));
                end
            end
        end

        function gTerm = finalizeRHSGravTerm(obj,lw)            
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'SPFlow'));
            nSubCells = nnz(obj.dofm.subCells(:,subInd));
            intFaces = any(obj.isIntFaces(:,subInd),2);
            neigh = obj.faces.faceNeighbors(intFaces,:);
            neigh = obj.dofm.ent2field('SPFlow',neigh);
            gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
                 -lw.*obj.rhsGrav],[nSubCells,1]);

        end

        function blk = blockJacobian(obj,varargin)
            fRow = varargin{1}; 
            fCol = varargin{2};
            dt = varargin{3};
            locRow = obj.dofm.field2block(fRow);
            locCol = obj.dofm.field2block(fCol);
            blk = obj.simParams.theta*obj.H(locRow,locCol) + obj.P(locRow,locCol)/dt;
        end

        function blk = blockRhs(obj, fld)
            if ~strcmp(obj.dofm.subPhysics(fld), 'SPFlow')
                % no contribution to non poro fields
                blk = 0;
            else
                dofs = obj.dofm.field2block(fld);
                blk = obj.rhs(dofs);
            end
        end

        function out = isLinear(obj)
            out = true;
        end


        function trans = getFaceTransmissibilities(obj,faceID)
            trans = obj.trans(faceID);
        end

        function computeTrans(obj)   % Inspired by MRST
            % Compute first the vector connecting each cell centroid to the
            % half-face
            r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
            c = [1, 2, 3, 1, 2, 3, 1, 2, 3];
            hf2Cell = repelem((1:obj.mesh.nCells)',diff(obj.faces.mapF2E));
            L = obj.faces.faceCentroid(obj.faces.faces2Elements(:,1),:) - obj.elements.cellCentroid(hf2Cell,:);
            sgn = 2*(hf2Cell == obj.faces.faceNeighbors(obj.faces.faces2Elements(:,1))) - 1;
            N = bsxfun(@times,sgn,obj.faces.faceNormal(obj.faces.faces2Elements(:,1),:));
            KMat = zeros(obj.mesh.nCellTag,9);
            for i=1:obj.mesh.nCellTag
                KMat(i,:) = obj.material.getMaterial(i).PorousRock.getPermVector();
            end
            hT = zeros(length(hf2Cell),1);
            for k=1:length(r)
                hT = hT + L(:,r(k)) .* KMat(obj.mesh.cellTag(hf2Cell),k) .* N(:,c(k));
            end
            hT = hT./sum(L.*L,2);
            %       mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
            %       hT = hT/mu;
            %
            obj.trans = 1 ./ accumarray(obj.faces.faces2Elements(:,1),1 ./ hT,[obj.faces.nFaces,1]);
        end

    end
end

