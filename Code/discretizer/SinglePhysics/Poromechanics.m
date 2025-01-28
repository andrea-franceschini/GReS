classdef Poromechanics < SinglePhysicSolver
    %POROMECHANICS
    properties
        nEntryKLoc      % Entries of local stiffness matrix
        fInt            % internal forces
    end

    methods (Access = public)
        function obj = Poromechanics(symmod,params,dofManager,grid,mat,data)
            obj@SinglePhysicSolver('Poromechanics',symmod,params,dofManager,grid,mat,data);
            obj.nEntryKLoc = (obj.mesh.nDim^2)*(obj.elements.nNodesElem).^2;
        end

        function state = computeMat(obj,varargin)
            % Compute Stiffness matrix for mechanical problem
            n = sum(~cellfun(@isempty,varargin));
            if n == 2
              state = varargin{1};
              dt = varargin{2};
            elseif n == 3
              state = varargin{1};
              dt = varargin{3};
            end
            subCells = obj.dofm.getFieldCells(obj.field);
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iiVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            jjVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            KVec = zeros(obj.nEntryKLoc*nSubCellsByType,1);
            nDof = obj.dofm.getNumDoF(obj.field);
            obj.fInt = zeros(nDof,1); % reset internal forces
            nComp = obj.dofm.getDoFperEnt(obj.field);
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
                dof = dofId(obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)),nComp);
                [jjLoc,iiLoc] = meshgrid(dof,dof);
                iiVec(l1+1:l1+s1) = iiLoc(:);
                jjVec(l1+1:l1+s1) = jjLoc(:);
                KVec(l1+1:l1+s1) = KLoc(:);
                l1 = l1 + s1;
                l2 = l2 + s2;
                obj.fInt(dof) = obj.fInt(dof)+fLoc;
            end
            % renumber indices according to active nodes
            [~,~,iiVec] = unique(iiVec);
            [~,~,jjVec] = unique(jjVec);
            % populate stiffness matrix
            obj.J = sparse(iiVec, jjVec, KVec, nDof, nDof);
            if obj.simParams.isTimeDependent
               obj.J = obj.simParams.theta*obj.J;
            end
        end

        function state = setState(obj,state)
           % add poromechanics fields to state structure
           ng = 1;
           if ~isempty(obj.GaussPts)
              ng = obj.GaussPts.nNode;
           end
           state.conv = struct('strain', [], 'stress', [], 'status', []);
           state.curr = struct('strain', [], 'stress', [], 'status', []);
           state.curr.stress = zeros([1, ng, 0, 0]*obj.elements.nCellsByType,6);
           state.curr.strain = zeros([1, ng, 0, 0]*obj.elements.nCellsByType,6);
           state.curr.status = zeros([1, ng, 0, 0]*obj.elements.nCellsByType,2);
           state.conv = state.curr;
           state.iniStress = zeros([1, ng, 0, 0]*obj.elements.nCellsByType,6);
           state.dispConv = zeros(obj.mesh.nDim*obj.mesh.nNodes,1);
           state.dispCurr = zeros(obj.mesh.nDim*obj.mesh.nNodes,1);
        end

        function state = updateState(obj,state,dSol)
           % Update state structure with last solution increment
           ents = obj.dofm.getActiveEnts(obj.field);
           state.dispCurr(ents) = state.dispCurr(ents) + dSol(getDoF(obj.dofm,obj.field));
           du = state.dispCurr - state.dispConv;
           % Update strain
           l = 0;
           for el=1:obj.mesh.nCells
              dof = getDoFID(obj.mesh,el);
              switch obj.mesh.cellVTKType(el)
                 case 10 % Tetra
                    N = getDerBasisF(obj.elements.tetra,el);
                    B = zeros(6,4*obj.mesh.nDim);
                    B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
                    state.curr.strain(l+1,:) = (B*du(dof))';
                    l = l + 1;
                 case 12 % Hexa
                    N = getDerBasisFAndDet(obj.elements.hexa,el,2);
                    B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                    B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                    state.curr.strain((l+1):(l+obj.GaussPts.nNode),:) = ...
                       reshape(pagemtimes(B,du(dof)),6,obj.GaussPts.nNode)';
                    l = l + obj.GaussPts.nNode;
              end
           end
        end

        function [avStress,avStrain] = finalizeState(obj,state)
           % compute cell average values of stress and strain
           avStress = zeros(obj.mesh.nCells,6);
           avStrain = zeros(obj.mesh.nCells,6);
           l = 0;
           for el=1:obj.mesh.nCells
              dof = getDoFID(obj.mesh,el);
              switch obj.mesh.cellVTKType(el)
                 case 10 % Tetra
                    N = getDerBasisF(obj.elements.tetra,el);
                    B = zeros(6,4*obj.mesh.nDim);
                    B(obj.elements.indB(1:36,2)) = N(obj.elements.indB(1:36,1));
                    dStrain = B*state.dispCurr(dof);
                    avStrain(el,:) = dStrain;
                    avStress(el,:) = state.conv.stress(l+1,:);
                    l = l + 1;
                 case 12 % Hexa
                    [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                    B = zeros(6,8*obj.mesh.nDim,obj.GaussPts.nNode);
                    B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                    avStress(el,:) = sum(diag(dJWeighed)* ...
                       state.curr.stress((l+1):(l+obj.GaussPts.nNode),:))/obj.elements.vol(el);
                    dStrain = pagemtimes(B,state.dispCurr(dof));
                    dStrain = dStrain.*reshape(dJWeighed,1,1,[]);
                    avStrain(el,:) = sum(dStrain,3)/obj.elements.vol(el);
                    l = l + obj.GaussPts.nNode;
              end
           end
        end


        function stateTmp = computeRhs(obj,varargin)
           stateTmp = varargin{1};
           %computing rhs for poromechanics field
           if isLinear(obj) % linear case
              % update elastic stress
              subCells = obj.dofm.getFieldCells(obj.field);
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
              ents = obj.dofm.getActiveEnts(obj.field);
              if obj.simParams.isTimeDependent
                 theta = obj.simParams.theta;
                 obj.rhs = obj.J*stateTmp.dispCurr(ents) + ...
                    (1/theta-1)*obj.J*stateTmp.dispConv(ents);
              else
                 obj.rhs = obj.J*stateTmp.dispCurr(ents);
              end
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

