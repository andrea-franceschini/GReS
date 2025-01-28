classdef SPFlow < SinglePhysicSolver
    %POROMECHANICS
    % Subclass of Discretizer
    % Implement Poromechanics methods to assemble the stiffness matrix and
    % the residual contribution

    properties
        trans
        isIntFaces
        rhsGrav         % gravity contribution to rhs
        H
        P
    end

    methods (Access = public)
        function obj = SPFlow(symmod,params,dofManager,grid,mat,data)
            obj@SinglePhysicSolver('SPFlow',symmod,params,dofManager,grid,mat,data);
            if obj.model.isFVTPFABased('Flow')
                obj.computeTrans;
                %get cells with active flow model
                flowCells = obj.dofm.getFieldCells(obj.field);
                % Find internal faces (i.e. shared by two cells with active
                % flow)
                obj.isIntFaces = all(ismember(obj.faces.faceNeighbors, flowCells), 2);
            end
            computeRHSGravTerm(obj);
        end

        function state = setState(obj,state)
           if isFEMBased(obj.model,'Flow')
              n = obj.mesh.nNodes;
           elseif isFVTPFABased(obj.model,'Flow')
              n = obj.mesh.nCells;
           end
           state.pressure = zeros(n,1);
        end

        function updateState(obj,state,dSol)
           ents = obj.dofm.getActiveEnts(obj.field);
           state.pressure(ents) = state.pressure(ents) + dSol(obj.dofm.getDoF(obj.field));
        end

        function fluidPot = finalizeState(obj,state)
           fluidPot = state.pressure;
           gamma = obj.material.getFluid().getFluidSpecWeight();
           if gamma > 0
              if isFEMBased(obj.model,'Flow')
                 fluidPot = fluidPot + gamma*obj.mesh.coordinates(:,3);
              elseif isFVTPFABased(obj.model,'Flow')
                 fluidPot = fluidPot + gamma*obj.elements.cellCentroid(:,3);
              end
           end
        end


        function state = computeMat(obj,varargin)
           state = varargin{1};
           dt = varargin{3};
            if obj.model.isFEMBased('Flow')
                computeMatFEM(obj);
            elseif obj.model.isFVTPFABased('Flow')
                mu = obj.material.getFluid().getDynViscosity();
                computeStiffMatFV(obj,1/mu);
                computeCapMatFV(obj);
            end
            if obj.simParams.isTimeDependent
               obj.J = obj.simParams.theta*obj.H + obj.P/dt;
            else
               obj.J = obj.H;
            end
        end

        function computeMatFEM(obj) 
            % dealing with input params
            subCells = obj.dofm.getFieldCells(obj.field);
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
                dof = dofId(obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)),1);
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                HVec(l1+1:l1+s1) = HLoc(:);
                PVec(l1+1:l1+s1) = PLoc(:);
                l1 = l1 + s1;
            end
            % renumber indices according to active nodes
            [~,~,iiVec] = unique(iiVec);
            [~,~,jjVec] = unique(jjVec);
            nDoF = obj.dofm.getNumDoF(obj.field);
            % Assemble H and P matrices defined as new fields of
            obj.H = sparse(iiVec, jjVec, HVec, nDoF, nDoF);
            obj.P = sparse(iiVec, jjVec, PVec, nDoF, nDoF);
        end


        function computeStiffMatFV(obj,lw)
            % Inspired by MRST
            subCells = obj.dofm.getFieldCells('Poromechanics');
            nSubCells = length(subCells); 
            %get pairs of faces that contribute to the subdomain
            neigh1 = obj.faces.faceNeighbors(obj.isIntFaces,1);
            neigh2 = obj.faces.faceNeighbors(obj.isIntFaces,2);
            % Transmissibility of internal faces
            tmpVec = lw.*obj.trans(obj.isIntFaces);
            % tmpVec = lw.*tmpVec;
            sumDiagTrans = accumarray([neigh1; neigh2], ...
                repmat(tmpVec,[2,1]),[nSubCells,1]);
            % Assemble H matrix
            nDoF = obj.dofm.getNumDoF(obj.field);
            [~,~,neigh1] = unique(neigh1); 
            [~,~,neigh2] = unique(neigh2); 
            obj.H = sparse([neigh1; neigh2; (1:nSubCells)'], ...
                [neigh2; neigh1; (1:nSubCells)'], ...
                [-tmpVec; -tmpVec; ...
                sumDiagTrans],nDoF,nDoF);
        end



        function computeCapMatFV(obj,varargin)
            subCells = obj.dofm.getFieldCells(obj.field);
            nSubCells = length(subCells);
            poroMat = zeros(nSubCells,1);
            alphaMat = zeros(nSubCells,1);
            beta = obj.material.getFluid().getFluidCompressibility();
            for m = 1:obj.mesh.nCellTag
                if ~ismember(m,obj.dofm.getFieldCellTags({obj.field,'Poromechanics'}))
                    % compute alpha only if there's no coupling in the
                    % subdomain
                    alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
                end
                poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
            end
            % (alpha+poro*beta)
            PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));
            PVal = PVal.*obj.elements.vol(subCells);
            nDoF = obj.dofm.getNumDoF(obj.field);
            [~,~,dof] = unique(subCells);
            obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
        end


        function stateTmp = computeRhs(obj,stateTmp,statek,dt)
            % Compute the residual of the flow problem
            lw = obj.material.getFluid().getDynViscosity();
            ents = obj.dofm.getActiveEnts(obj.field);
            if obj.simParams.isTimeDependent
               obj.rhs = obj.H*stateTmp.pressure(ents);
            else
               theta = obj.simParams.theta;
               rhsStiff = theta*obj.H*stateTmp.pressure(ents) + (1-theta)*obj.H*statek.pressure(ents);
               rhsCap = (obj.P/dt)*(stateTmp.pressure(ents) - statek.pressure(ents));
               obj.rhs = rhsStiff + rhsCap;
            end
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
            rhsTmp = zeros(obj.dofm.getNumDoF(obj.field),1);
            gamma = obj.material.getFluid().getFluidSpecWeight();
            if gamma > 0
                subCells = obj.dofm.getFieldCells(obj.field);
                if isFEMBased(obj.model,'Flow')
                    for el = subCells'
                        % Get the material permeability
                        permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
                        %             permMat = permMat/mu;
                        switch obj.mesh.cellVTKType(el)
                            case 10 % Tetrahedra
                                N = obj.elements.tetra.getDerBasisF(el);
                                rhsLoc = (N'*permMat(:,3))*obj.elements.vol(el)*gamma;
                            case 12 % Hexa
                                [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
                                fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
                                fs = fs.*reshape(dJWeighed,1,1,[]);
                                rhsLoc = sum(fs,3)*gamma;
                        end
                        %
                        entsId = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
                        rhsTmp = rhsTmp(entsId) + rhsLoc; 
                    end
                    obj.rhsGrav = rhsTmp(obj.dofm.getActiveEnts(obj.field));
                elseif isFVTPFABased(obj.model,'Flow')
                    neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
                    zVec = obj.elements.cellCentroid(:,3);
                    zNeigh = zVec(neigh);
                    obj.rhsGrav = gamma*obj.trans(obj.isIntFaces).*(zNeigh(:,1) - zNeigh(:,2));
                end
            end
            % remove inactive components of rhs vector 
        end

        function gTerm = finalizeRHSGravTerm(obj,lw) 
            nCells = obj.dofm.getNumDoF(obj.field);
            neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
            gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
                 -lw.*obj.rhsGrav],[nCells,1]);
            gTerm = gTerm(obj.dofm.getActiveEnts(obj.field));
        end

        function applyDirBC(obj,~,ents,vals,state)
           % apply Dirichlet BCs
           % ents: id of constrained faces without any dof mapping applied
           if isFVTPFABased(obj.model,'Flow')
              % BCs imposition for finite volumes
              neigh = sum(obj.faces.faceNeighbors(ents,:),2);
              % deal with cells having > 1 constrained faces
              [neigh,~,ind] = unique(neigh);    
              gamma = obj.material.getFluid().getFluidSpecWeight();
              mu = obj.material.getFluid().getDynViscosity();
              tr = obj.getFaceTransmissibilities(ents);
              q = 1/mu*tr.*((state.pressure(neigh) - vals)...
                 + gamma*(obj.elements.cellCentroid(neigh,3) - obj.faces.faceCentroid(ents,3)));
              rhsVals = accumarray(ind,q);
              JVals = 1/mu*tr;
              dofs = obj.dofm.getLocalDoF(neigh,obj.field); % get dof numbering 
              nDoF = obj.dofm.getNumDoF(obj.field);
              obj.J(nDoF*(dofs-1) + nDoF) = obj.J(nDoF*(dofs-1) + nDoF) + JVals;
              obj.rhs(dofs) = obj.rhs(dofs) + rhsVals;
           else
              applyDirBC@SinglePhysicSolver(obj,ents)
           end
        end

        function applyNeuBC(obj,ents,vals)
           if isFVTPFABased(obj.model,'Flow')
           area = vecnorm(obj.faces.faceNormal(ents,:),2,2).*vals;
           rhsVals = accumarray(ind, area);
           end
           applyNeuBC@SinglePhysicSolver(obj,ents,rhsVals)
        end

        function applyVolumeForceBC(obj,dofs,vals)
           % map local dof numbering to global entitities numbering
           entID = obj.dofm.getFieldDoF(dofs,obj.field);
           vols = obj.elements.vol(entID);
           rhsVals = vals.*vols;
           obj.rhs(dofs) = obj.rhs(dofs)- rhsVals;
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

