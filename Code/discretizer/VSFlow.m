classdef VSFlow < handle
    % Variably Saturated flow
    % Subclass of Discretizer
    % Implements Richards equations for unsaturated flow in vadose region
    % most of the method are directly derived from SPFlow class, but the
    % class is different to preserve code modularity

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
        lwkpt           % mobility
        H               % stiffness matrix
        P               % capacity matrix
        JNewt = []      % newton jacobian contribution
        rhs             % residual
        rhsGrav         % gravity contribute to rhs
        upElem          % upstream elements array for each face
    end

    methods (Access = public)
        function obj = VSFlow(symmod,params,dofManager,grid,mat,data)
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
                    if any(strcmp(obj.dofm.subDomains(i).physics,"VSFlow"))
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
                obj.upElem = zeros(nnz(obj.isIntFaces),1);
            end
            computeRHSGravTerm(obj);
        end



        function computeMat(obj,varargin)
            stateTmp = varargin{1};
            statek = varargin{2};
            dt = varargin{3};
            pkpt = obj.simParams.theta*stateTmp.pressure + ...
                (1 - obj.simParams.theta)*statek.pressure;
            [Swkpt,dSwkpt,obj.lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt);
            computeStiffMatFV(obj,obj.lwkpt);
            computeCapMatFV(obj,Swkpt,dSwkpt);
            if isNewtonNLSolver(obj.simParams)
                computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,dlwkpt)
            end
        end


        function computeStiffMatFV(obj,lw)
            % Inspired by MRST
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCells = length(subCells); 
            %get pairs of faces that contribute to the subdomain
            intFaces = any(obj.isIntFaces(:,subInd),2);
            neigh1 = obj.faces.faceNeighbors(intFaces,1);
            neigh2 = obj.faces.faceNeighbors(intFaces,2);
            neigh1dof = obj.dofm.ent2field('VSFlow',neigh1);
            neigh2dof = obj.dofm.ent2field('VSFlow',neigh2);
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
            % the model is uncoupled
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
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
            % compute storage using updated saturation
            PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));     
            PVal = PVal.*varargin{1} + poroMat(obj.mesh.cellTag).*varargin{2};
            PVal = PVal.*obj.elements.vol(subCells);
            dof = obj.dofm.ent2field('VSFlow', subCells);
            obj.P = sparse(dof,dof,PVal);
        end


        function computeRhs(obj,stateTmp,statek,dt)
            % Compute the residual of the flow problem
            theta = obj.simParams.theta;
            ents = obj.dofm.field2ent('VSFlow');
            rhsStiff = theta*obj.H*stateTmp.pressure(ents) + (1-theta)*obj.H*statek.pressure(ents);
            rhsCap = (obj.P/dt)*(stateTmp.pressure(ents) - statek.pressure(ents));
            obj.rhs = rhsStiff + rhsCap;
            gamma = obj.material.getFluid().getFluidSpecWeight();
            %adding gravity rhs contribute
            if gamma > 0
                if isFEMBased(obj.model,'Flow')
                    obj.rhs = obj.rhs + obj.rhsGrav;
                elseif isFVTPFABased(obj.model,'Flow')
                    obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,obj.lwkpt);
                end
            end
        end

        function computeRHSGravTerm(obj)
            % Compute the gravity contribution
            % Get the fluid specific weight and viscosity'
            rhsTmp = zeros(sum(obj.dofm.numDof),1);
            gamma = obj.material.getFluid().getFluidSpecWeight();
            if gamma > 0
                subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
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
                        dof = obj.dofm.ent2field('VSFlow',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                        rhsTmp(dof) = rhsTmp(dof) + fLoc; 
                    end
                    obj.rhsGrav = rhsTmp(obj.dofm.getDoF('VSFlow'));
                elseif isFVTPFABased(obj.model,'Flow')
                    intFaces = any(obj.isIntFaces(:,subInd),2);
                    neigh = obj.faces.faceNeighbors(intFaces,:);
                    zVec = obj.elements.cellCentroid(:,3);
                    zNeigh = zVec(neigh);
                    obj.rhsGrav = gamma*obj.trans(intFaces).*(zNeigh(:,1) - zNeigh(:,2));
                end
            end
        end

        function gTerm = finalizeRHSGravTerm(obj,lw)            
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
            nSubCells = nnz(obj.dofm.subCells(:,subInd));
            intFaces = any(obj.isIntFaces(:,subInd),2);
            neigh = obj.faces.faceNeighbors(intFaces,:);
            neigh = obj.dofm.ent2field('VSFlow',neigh(:));
            gTerm = accumarray(neigh,[lw.*obj.rhsGrav; ...
                 -lw.*obj.rhsGrav],[nSubCells,1]);

        end

        function blk = blockJacobian(obj,varargin)
            fRow = varargin{1}; 
            fCol = varargin{2};
            dt = varargin{3};
            locRow = obj.dofm.field2block(fRow);
            locCol = obj.dofm.field2block(fCol);
            if isNewtonNLSolver(obj.simParams)
                blk = obj.simParams.theta*obj.H(locRow,locCol) + obj.P(locRow,locCol)/dt + ...
                    obj.JNewt(locRow,locCol);
            else
                blk = obj.simParams.theta*obj.H(locRow,locCol) + obj.P(locRow,locCol)/dt;

            end
        end

        function blk = blockRhs(obj, fld)
            if ~strcmp(obj.dofm.subPhysics(fld), 'VSFlow')
                % no contribution to non poro fields
                blk = 0;
            else
                dofs = obj.dofm.field2block(fld);
                blk = obj.rhs(dofs);
            end
        end

        function out = isLinear(obj)
            out = false;
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

        function [Swkpt,dSwkpt,lwkpt,dlwkpt] = computeUpElemAndProperties(obj,pkpt)
            % compute upstream elements for each face
            % interpolate effective saturation and relative permeability 
            % and first derivatives
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
            intFaces = any(obj.isIntFaces(:,subInd),2);
            neigh = obj.faces.faceNeighbors(intFaces,:);
            gamma = obj.material.getFluid().getFluidSpecWeight();
            if gamma > 0
                zVec = obj.elements.cellCentroid(:,3);
                zNeigh = zVec(neigh);
                lElemIsUp = pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)) >= 0;
            else
                lElemIsUp = pkpt(neigh(:,1)) >= pkpt(neigh(:,2));
            end
            obj.upElem(lElemIsUp) = neigh(lElemIsUp,1);
            obj.upElem(~lElemIsUp) = neigh(~lElemIsUp,2);
            [Swkpt,dSwkpt] = obj.material.computeSwAnddSw(obj.mesh,pkpt);
            dSwkpt = - dSwkpt;
            [lwkpt,dlwkpt] = obj.material.computeLwAnddLw(obj.mesh,obj.upElem,pkpt);
            dlwkpt = - dlwkpt;
        end

        function computeNewtPartOfJacobian(obj,dt,statek,stateTmp,pkpt,dSwkpt,dlwkpt)
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'VSFlow'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCells = length(subCells);
            intFaces = any(obj.isIntFaces(:,subInd),2);
            % compute matrices J1 and J2 (gathering non linear terms)
            neigh = obj.faces.faceNeighbors(intFaces,:);
            zVec = obj.elements.cellCentroid(:,3);
            zNeigh = zVec(neigh);
            gamma = obj.material.getFluid().getFluidSpecWeight();
            tmpVec1 = (dlwkpt.*obj.trans(intFaces)).*(pkpt(neigh(:,1)) - pkpt(neigh(:,2)) + gamma*(zNeigh(:,1) - zNeigh(:,2)));
            %
            poroMat = zeros(obj.mesh.nCellTag,1);
            alphaMat = zeros(obj.mesh.nCellTag,1);
            beta = obj.material.getFluid().getFluidCompressibility();
            for m = 1:obj.mesh.nCellTag
                poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
                alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
            end
            tmpVec2 = alphaMat(obj.mesh.cellTag) + beta*poroMat(obj.mesh.cellTag);
            tmpVec2 = 1/dt*((tmpVec2(subCells).*dSwkpt(subCells)).*(stateTmp.pressure(subCells) - statek.pressure(subCells))).*obj.elements.vol(subCells);
            neigh1dof = obj.dofm.ent2field('VSFlow',neigh(:,1));
            neigh2dof = obj.dofm.ent2field('VSFlow',neigh(:,2));
            upElemdof = obj.dofm.ent2field('VSFlow',obj.upElem);
            obj.JNewt = sparse([neigh1dof; neigh2dof; (1:nSubCells)'], ...
                [repmat(upElemdof,[2,1]);  (1:nSubCells)'], ...
                [tmpVec1; -tmpVec1; tmpVec2],nSubCells,nSubCells);
            obj.JNewt = obj.simParams.theta*obj.JNewt;
        end
    end
end

