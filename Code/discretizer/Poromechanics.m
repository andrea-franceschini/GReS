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
        fInt            % internal forces
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
            n = sum(~cellfun(@isempty,varargin));
            if n == 2
              state = varargin{1};
              dt = varargin{2};
            elseif n == 3
              state = varargin{1};
              dt = varargin{3};
            end
            subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'Poro'));
            [subCells, ~] = find(obj.dofm.subCells(:,subInd));
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iiVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            jjVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            KVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            nDof = obj.dofm.getNumbDof('Poro');
            obj.fInt = zeros(nDof,1); % reset internal forces
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
                        if ~isLinear(obj)
                           state.curr.stress(l2+1,:) = sigma;
                        end
                        KLoc = B'*D*B*vol;
                        s1 = obj.nEntryKLoc(1);
                        sz = sigma - state.iniStress(l2+1,:);
                        fLoc = (B')*sz'*vol;
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
                        fTmp = pagemtimes(B,'ctranspose',sz,'none');
                        fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
                        fLoc = sum(fTmp,3);
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
                obj.fInt(dof) = obj.fInt(dof)+fLoc;
            end
            % % populate stiffness matrix
            obj.K = sparse(iiVec, jjVec, KVec);
        end


        function computeRhs(obj,varargin)
           stateTmp = varargin{1};
           %computing rhs for poromechanics field
           if isLinear(obj) % linear case
              % update elastic stress
              subInd = obj.dofm.subList(ismember(obj.dofm.subPhysics, 'Poro'));
              [subCells, ~] = find(obj.dofm.subCells(:,subInd));
              l1 = 0;
              for el=subCells'
                 D = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw.Dmat;
                 % Get the right material stiffness for each element
                 switch obj.mesh.cellVTKType(el)
                    case 10 % Tetrahedra
                       stateTmp.curr.stress(l1+1,:) = stateTmp.curr.stress(l1+1,:)+stateTmp.curr.strain(l1+1,:)*D;
                       s1 = 1;
                    case 12 % Hexahedra
                       stateTmp.curr.stress((l1+1):(l1+obj.GaussPts.nNode),:) = ...
                          stateTmp.curr.stress((l1+1):(l1+obj.GaussPts.nNode),:)+...
                          stateTmp.curr.strain((l1+1):(l1+obj.GaussPts.nNode),:)*D;
                       s1 = obj.GaussPts.nNode;
                 end
                 l1 = l1+s1;
              end
                 ents = obj.dofm.field2ent('Poro');
                 theta = obj.simParams.theta;
                 obj.rhs = theta*obj.K*stateTmp.dispCurr(ents) + ...
                  (1-theta)*obj.K*stateTmp.dispConv(ents);
                
            else % non linear case: rhs computed with internal forces (B^T*sigma)
                obj.rhs = obj.fInt; % provisional assuming theta = 1;
           end
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
            out = true;
              for i = 1:obj.mesh.nCellTag
                out = isa(obj.material.getMaterial(i).ConstLaw,"Elastic");
                if ~out
                  return;
                end
              end
        end
    end
end
