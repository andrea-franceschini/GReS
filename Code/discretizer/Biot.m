classdef Biot < handle
    % Biot model subclass
    % Coupled poromechanics with:
    % SinglePhase flow 
    % Variably Saturated flow: still to be updated

    properties
        model
        simParams
        dofm
        mesh
        elements
        faces
        material
        GaussPts
        flowStr
        Q               % stiffness matrix
        rhsPoro         % residual contribution to Poro problem
        rhsFlow         % residual contribution to Flow problem
    end

    methods (Access = public)
        function obj = Biot(symmod,params,dofManager,grid,mat,data)
            %POROMECHANICS Construct an instance of this class
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

            if isSinglePhaseFlow(obj.model)
                obj.flowStr = 'SPFlow';
            elseif isVariabSatFlow(obj.model)
                obj.flowStr = 'VSFlow';
            end
            %
        end

        function computeMat(obj,varargin)
            if isFEMBased(obj.model, 'Flow')
                computeMatFEM_FEM(obj);
            elseif isFVTPFABased(obj.model,'Flow')
                computeMatFEM_FV(obj);
            end
        end


        function computeMatFEM_FEM(obj)
            % get domain ID with active Poro and Flow field
            nSub = length(obj.dofm.subDomains);
            subInd = zeros(nSub,1);
            for i = 1:length(obj.dofm.subDomains)
                if all(ismember(["Poro"; obj.flowStr], obj.dofm.subDomains(i).physics))
                    subInd(i) = i;
                end
            end
            subInd = sort(nonzeros(subInd));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            Qvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            %
            l1 = 0;
            if nSubCellsByType(2) > 0
                N1 = obj.elements.hexa.getBasisFinGPoints();
            end
            for el=subCells'
                % Get the right material stiffness for each element
                biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
                switch obj.mesh.cellVTKType(el)
                    case 10 %Tetrahedrons, direct integration
                        vol = findVolume(obj.elements.tetra,el);
                        der = getDerBasisF(obj.elements.tetra,el);
                        Qloc = biot*0.25*repelem(der(:),1,4)*vol;
                        s1 = obj.elements.nNodesElem(1).^2*obj.mesh.nDim;
                    case 12 %Hexahedrons, Gauss integration
                        [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                        nG = obj.GaussPts.nNode;
                        iN = zeros(6,8,nG); %matrix product i*N
                        B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                        B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                        iN(1,:,:) = reshape(N1',1,8,nG);
                        iN(2:3,:,:) = repmat(iN(1,:,:),2,1,1);
                        Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
                        Qs = Qs.*reshape(dJWeighed,1,1,[]);
                        Qloc = sum(Qs,3);
                        clear Qs;
                        s1 = obj.elements.nNodesElem(2)^2*obj.mesh.nDim;
                end
                %
                %assembly Coupling Matrix
                dofrow = obj.dofm.ent2field('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                dofcol = obj.dofm.ent2field(obj.flowStr,obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                iivec(l1+1:l1+s1) = iiloc(:);
                jjvec(l1+1:l1+s1) = jjloc(:);
                Qvec(l1+1:l1+s1) = Qloc(:);
                l1 = l1+s1;
            end
            nDofPoro = sum(obj.dofm.numDof(strcmp(obj.dofm.subPhysics,'Poro')));
            nDofFlow = sum(obj.dofm.numDof(strcmp(obj.dofm.subPhysics,obj.flowStr)));
            obj.Q = sparse(iivec, jjvec, Qvec, nDofPoro, nDofFlow);
        end

        function computeMatFEM_FV(obj)
            % get domain ID with active Poro and Flow field
            nSub = length(obj.dofm.subDomains);
            subInd = zeros(nSub,1);
            for i = 1:length(obj.dofm.subDomains)
                if all(ismember(["Poro"; obj.flowStr], obj.dofm.subDomains(i).physics))
                    subInd(i) = i;
                end
            end
            subInd = nonzeros(subInd);
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            Qvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            %
            l1 = 0;
            if nSubCellsByType(2) > 0
                N1 = obj.elements.hexa.getBasisFinGPoints();
            end
            for el=subCells'
                % Get the right material stiffness for each element
                biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
                switch obj.mesh.cellVTKType(el)
                    case 10 %Tetrahedrons, direct integration
                        error('TPFA-FV not implemented for Tetrahedrons')
                        % vol = findVolume(obj.elements.tetra,el);
                        % der = getDerBasisF(obj.elements.tetra,el);
                        % Qloc = biot*0.25*repelem(der(:),1,4)*vol;
                        % s1 = obj.elements.nNodesElem(1).^2*obj.mesh.nDim;
                    case 12 %Hexahedrons, Gauss integration
                        [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                        nG = obj.GaussPts.nNode;
                        %iN = zeros(6,8,nG); %matrix product i*N
                        B = zeros(6,8*obj.mesh.nDim,nG);
                        B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                        tmp = [1;1;1;0;0;0];
                        kron = repmat(tmp,1,1,8);
                        % iN(1,:,:) = reshape(N1',1,8,nG);
                        % iN(2:3,:,:) = repmat(iN(1,:,:),2,1,1);
                        Qs = biot*pagemtimes(B,'ctranspose',kron,'none');
                        Qs = Qs.*reshape(dJWeighed,1,1,[]);
                        Qloc = sum(Qs,3);
                        clear Qs;
                        s1 = obj.elements.nNodesElem(2)*obj.mesh.nDim;
                end
                %
                %assembly Coupling Matrix
                dofrow = obj.dofm.ent2field('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                dofcol = obj.dofm.ent2field(obj.flowStr,el);
                [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                iivec(l1+1:l1+s1) = iiloc(:);
                jjvec(l1+1:l1+s1) = jjloc(:);
                Qvec(l1+1:l1+s1) = Qloc(:);
                l1 = l1+s1;
            end
            nDofPoro = sum(obj.dofm.numDof(strcmp(obj.dofm.subPhysics,'Poro')));
            nDofFlow = sum(obj.dofm.numDof(strcmp(obj.dofm.subPhysics,obj.flowStr)));
            obj.Q = sparse(iivec, jjvec, Qvec, nDofPoro, nDofFlow);
        end

        function computeRhs(obj,stateTmp,statek,dt)
            % compute Biot rhs contribute
            theta = obj.simParams.theta;
            % select active coefficients of matrices and solution vectors
            entsPoro = obj.dofm.field2ent('Poro');
            entsFlow = obj.dofm.field2ent(obj.flowStr);
            % maybe im taking into account displacement nodes that shuold
            % not be coupled!
            %
            obj.rhsPoro = -theta*obj.Q*stateTmp.pressure(entsFlow) - ...
                (1-theta)*(obj.Q*statek.pressure(entsFlow));
            obj.rhsFlow = (obj.Q/dt)'*(stateTmp.dispCurr(entsPoro) - stateTmp.dispConv(entsPoro));
        end


        function blk = blockJacobian(obj,varargin)
            fRow = varargin{1}; 
            fCol = varargin{2};
            dt = varargin{3};
            locRow = obj.dofm.field2block(fRow);
            locCol = obj.dofm.field2block(fCol);
            if strcmp(obj.dofm.subPhysics(fRow), 'Poro')
                blk = -obj.simParams.theta*obj.Q(locRow,locCol);
            elseif strcmp(obj.dofm.subPhysics(fRow), obj.flowStr)
                blk = (obj.Q(locCol,locRow))'/dt;
            end
        end



        function blk = blockRhs(obj, fld)
            dofs = obj.dofm.field2block(fld);
            if strcmp(obj.dofm.subPhysics(fld), 'Poro')
                blk = obj.rhsPoro(dofs);
            elseif strcmp(obj.dofm.subPhysics(fld), obj.flowStr)
                blk = obj.rhsFlow(dofs);
            end
        end

        function out = isLinear(obj)
            out = true;
        end


    end

end


