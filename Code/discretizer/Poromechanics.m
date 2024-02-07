classdef Poromechanics < handle
    %POROMECHANICS
    % Subclass of Discretizer
    % Implement Poromechanics methods to assemble the stiffness matrix and
    % the residual contribution

    properties
        nEntryKLoc
        model
        simParams
        dofm
        mesh
        elements
        faces
        material
        GaussPts
        K               % stiffness matrix
        rhs             % residual
    end

    methods (Access = public)
        function obj = Poromechanics(symmod,params,dofManager,grid,mat,data)
            %POROMECHANICS Construct an instance of this class
            obj.model = symmod;
            obj.simParams = params;
            obj.dofm = dofManager;
            obj.mesh = grid.topology;
            obj.elements = grid.cells;
            obj.faces = grid.faces;
            obj.material = mat;
            %       obj.probType = pType;
            %       obj.bound = bc;
            %       obj.BCName = BCName;
            %       obj.state = stat;
            if ~isempty(data)
                obj.GaussPts = data{1};
            end
            %
            obj.nEntryKLoc = (obj.mesh.nDim^2)*(obj.elements.nNodesElem).^2;
        end

        function computeMat(obj,varargin)
            % Compute Stiffness matrix for mechanical problem 
            state = varargin{1};
            dt = varargin{2};
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'Poro'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iiVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            jjVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            KVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            %
            l1 = 0;
            l2 = 0;
            for el=subCells'
                % Get the right material stiffness for each element
                switch obj.mesh.cellVTKType(el)
                    case 10 % Tetrahedra
                        N = getDerBasisF(obj.elements.tetra,el);
                        vol = findVolume(obj.elements.tetra,el);
                        B = zeros(6,4*obj.mesh.nDim);
                        B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
                        [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                            state.conv.stress(l2+1,:), ...
                            state.curr.strain(l2+1,:), ...
                            dt,state.conv.status(l2+1,:), el, state.t);
                        state.curr.status(l2+1,:) = status;
                        state.curr.stress(l2+1,:) = sigma;
                        % D = obj.preP.getStiffMatrix(el,state.stress(l2+1,3)+state.iniStress(l2+1,3));
                        KLoc = B'*D*B*vol;
                        s1 = obj.nEntryKLoc(1);
                        %
                        %             if obj.flCompRHS
                        sz = sigma - state.iniStress(l2+1,:);
                        %fLoc = (B')*sz'*vol;
                        s2 = 1;
                        %             end
                    case 12 % Hexahedra
                        [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                        B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                        B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                        [D, sigma, status] = obj.material.updateMaterial(obj.mesh.cellTag(el), ...
                            state.conv.stress(l2+1:l2+obj.GaussPts.nNode,:), ...
                            state.curr.strain(l2+1:l2+obj.GaussPts.nNode,:), ...
                            dt,state.conv.status(l2+1:l2+obj.GaussPts.nNode,:), el, state.t);
                        state.curr.status(l2+1:l2+obj.GaussPts.nNode,:) = status;
                        state.curr.stress((l2+1):(l2+obj.GaussPts.nNode),:) = sigma;
                        Ks = pagemtimes(pagemtimes(B,'ctranspose',D,'none'),B);
                        Ks = Ks.*reshape(dJWeighed,1,1,[]);
                        KLoc = sum(Ks,3);
                        clear Ks;
                        s1 = obj.nEntryKLoc(2);
                        sz = sigma - state.iniStress(l2+1:l2+obj.GaussPts.nNode,:);
                        sz = reshape(sz',6,1,obj.GaussPts.nNode);
                        % fTmp = pagemtimes(B,'ctranspose',sz,'none');
                        % fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
                        % fLoc = sum(fTmp,3);
                        s2 = obj.GaussPts.nNode;
                end
                % get global DoF
                dof = obj.dofm.ent2field('Poro',obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                KVec(l1+1:l1+s1) = KLoc(:);
                l1 = l1 + s1;
                l2 = l2 + s2;
            end
            % populate stiffness matrix
            obj.K = sparse(iiVec, jjVec, KVec);
        end


        function computeRhs(obj,varargin)
            stateTmp = varargin{1};
            % computing rhs for poromechanics field
            ents = obj.dofm.field2ent('Poro');
            theta = obj.simParams.theta;
            obj.rhs = theta*obj.K*stateTmp.dispCurr(ents) + ...
                (1-theta)*obj.K*stateTmp.dispConv(ents);
        end

        function blk = blockJacobian(obj,varargin)
            fRow = varargin{1};
            fCol = varargin{2};
            locRow = obj.dofm.field2block(fRow);
            locCol = obj.dofm.field2block(fCol);
            blk = obj.simParams.theta*obj.K(locRow,locCol);
        end

        function blk = blockRhs(obj, fld)
            if ~strcmp(obj.dofm.subPhysics(fld), 'Poro')
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


    end
end

