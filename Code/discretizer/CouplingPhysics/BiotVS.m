classdef BiotVS < handle
    % WORK IN PROGRESS

    properties
        model
        simParams
        dofm
        mesh
        elements
        faces
        material
        isIntFaces
        upElem
        GaussPts
        flowStr
        Q1              % coupling matrix (both in poro and flow eqs)
        Q2              % additional coupling term in Poro equation
        rhsPoro         % residual contribution to Poro problem
        rhsFlow         % residual contribution to Flow problem
    end

    methods (Access = public)
        function obj = BiotVS(symmod,params,dofManager,grid,mat,data)
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

        function computeMat(obj,varargin)
            stateTmp = varargin{1};
            statek = varargin{2};
            pkpt = obj.simParams.theta*stateTmp.pressure + ...
                (1 - obj.simParams.theta)*statek.pressure;
            [Swkpt,dSwkpt,~,~] = computeUpElemAndProperties(obj,pkpt);
            % call method according to the discretization technique chosen
            if isFEMBased(obj.model, 'Flow')
                error('FEM is not available for VSFlow')
            elseif isFVTPFABased(obj.model,'Flow')
                computeMatFV(obj, Swkpt, dSwkpt, pkpt);
            end
        end


        function computeMatFV(obj, Swkpt, dSwkpt, pkpt)
            % get domain ID with active Poro and Flow field
            nSub = length(obj.dofm.subDomains);
            subInd = zeros(nSub,1);
            for i = 1:length(obj.dofm.subDomains)
                if all(ismember(["Poro"; "VSFlow"], obj.dofm.subDomains(i).physics))
                    subInd(i) = i;
                end
            end
            subInd = nonzeros(subInd);
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iivec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            jjvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            Qvec1 = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            Qvec2 = Qvec1;
            %
            l1 = 0;
            for el=subCells'
                % Get the right material stiffness for each element
                biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
                switch obj.mesh.cellVTKType(el)
                    case 10 %Tetrahedrons, direct integration
                        error('TPFA-FV not implemented for Tetrahedrons')
                    case 12 %Hexahedrons, Gauss integration
                        [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                        nG = obj.GaussPts.nNode;
                        B = zeros(6,8*obj.mesh.nDim,nG);
                        B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                        tmp = [1;1;1;0;0;0];
                        kron = repmat(tmp,1,1,8); % Kronecker delta
                        % Q1 local contribute
                        Qs1 = biot*Swkpt(el)*pagemtimes(B,'ctranspose',kron,'none');
                        Qs1 = Qs1.*reshape(dJWeighed,1,1,[]);
                        Qloc1 = sum(Qs1,3);
                        clear Qs1;
                        % Q2 local contribute
                        Qs2 = biot*dSwkpt(el)*pkpt(el)*pagemtimes(B,'ctranspose',kron,'none');
                        Qs2 = Qs2.*reshape(dJWeighed,1,1,[]);
                        Qloc2 = sum(Qs2,3);
                        clear Qs2;
                        s1 = obj.elements.nNodesElem(2)*obj.mesh.nDim;
                end
                %
                %assembly Coupling Matrix
                dofrow = obj.dofm.ent2field('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                dofcol = obj.dofm.ent2field('VSFlow',el);
                [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                iivec(l1+1:l1+s1) = iiloc(:);
                jjvec(l1+1:l1+s1) = jjloc(:);
                Qvec1(l1+1:l1+s1) = Qloc1(:);
                Qvec2(l1+1:l1+s1) = Qloc2(:);
                l1 = l1+s1;
            end
            nDofPoro = sum(obj.dofm.numDof(strcmp(obj.dofm.subPhysics,'Poro')));
            nDofFlow = sum(obj.dofm.numDof(strcmp(obj.dofm.subPhysics,'VSFlow')));
            obj.Q1 = sparse(iivec, jjvec, Qvec1, nDofPoro, nDofFlow);
            obj.Q2 = sparse(iivec, jjvec, Qvec2, nDofPoro, nDofFlow);
        end

        function computeRhs(obj,stateTmp,statek,dt)
            % compute Biot rhs contribute
            theta = obj.simParams.theta;
            % select active coefficients of matrices and solution vectors
            entsPoro = obj.dofm.field2ent('Poro');
            entsFlow = obj.dofm.field2ent('VSFlow');
            % maybe im taking into account displacement nodes that shuold
            % not be coupled!
            %
            obj.rhsPoro = -theta*(obj.Q1+obj.Q2)*stateTmp.pressure(entsFlow) - ...
                (1-theta)*(obj.Q1+obj.Q2)*statek.pressure(entsFlow);
            obj.rhsFlow = (obj.Q1/dt)'*(stateTmp.dispCurr(entsPoro) - stateTmp.dispConv(entsPoro));
        end


        function blk = blockJacobian(obj,varargin)
            fRow = varargin{1}; 
            fCol = varargin{2};
            dt = varargin{3};
            locRow = obj.dofm.field2block(fRow);
            locCol = obj.dofm.field2block(fCol);
            if strcmp(obj.dofm.subPhysics(fRow), 'Poro')
                blk = -obj.simParams.theta*(obj.Q1(locRow,locCol)+obj.Q2(locRow,locCol));
            elseif strcmp(obj.dofm.subPhysics(fRow), 'VSFlow')
                blk = (obj.Q1(locCol,locRow))'/dt;
            end
        end



        function blk = blockRhs(obj, fld)
            dofs = obj.dofm.field2block(fld);
            if strcmp(obj.dofm.subPhysics(fld), 'Poro')
                blk = obj.rhsPoro(dofs);
            elseif strcmp(obj.dofm.subPhysics(fld), 'VSFlow')
                blk = obj.rhsFlow(dofs);
            end
        end

        function out = isLinear(obj)
            out = true;
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
            [Swkpt,dSwkpt,d2Swkpt] = obj.material.computeSwAnddSw(obj.mesh,pkpt);
            dSwkpt = - dSwkpt;
            [lwkpt,dlwkpt] = obj.material.computeLwAnddLw(obj.mesh,obj.upElem,pkpt);
            dlwkpt = - dlwkpt;
         end
    end

end


