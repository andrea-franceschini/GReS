classdef Biot < CouplingPhysics
    % Biot model subclass
    % Coupled poromechanics with Flow models
    % field1: Poromechanics
    % field2: SPFlow

    properties
        Q
        flowScheme   % Discretization scheme used for Flow
    end

    methods (Access = public)
        function obj = Biot(symmod,params,dofManager,grid,mat,data)
            obj@CouplingPhysics('Poromechanics','SPFlow',symmod,params,dofManager,grid,mat,data);
            if isSinglePhaseFlow(obj.model)
                obj.flowScheme = 'SPFlow';
            elseif isVariabSatFlow(obj.model)
                obj.flowScheme = 'VSFlow';
            end
            %
        end

        function state = computeMat(obj,varargin)
           state = varargin{1};
           dt = varargin{3};
           % call method according to the discretization technique chosen
           if isempty(obj.J{1}) || ~isLinear(obj)
              if isFEMBased(obj.model, 'Flow')
                 computeMatFEM_FEM(obj);
              elseif isFVTPFABased(obj.model,'Flow')
                 computeMatFEM_FV(obj);
              end
           end
           % J{1}: momentum balance   J{2}: mass balance equations
           obj.J{1} = -obj.simParams.theta*obj.Q;
           obj.J{2} = obj.Q'/dt;
        end


        function computeMatFEM_FEM(obj)
            % get domain ID with active Poro and Flow field
            subCells = obj.dofm.getFieldCells(obj.fields);
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iiVec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            jjVec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            Qvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem.^2)*nSubCellsByType,1);
            nDoF1 = obj.dofm.getNumDoF(obj.fields{1});
            nDoF2 = obj.dofm.getNumDoF(obj.fields{2});
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
                %assembly Coupling Matrix
                dofrow = getLocalDoF(obj.dofm,obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)),obj.fields{1});
                dofcol = getLocalDoF(obj.dofm,obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)),obj.fields{2});
                [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                iiVec(l1+1:l1+s1) = iiloc(:);
                jjVec(l1+1:l1+s1) = jjloc(:);
                Qvec(l1+1:l1+s1) = Qloc(:);
                l1 = l1+s1;
            end
            obj.Q = sparse(iiVec, jjVec, Qvec, nDoF1, nDoF2);
        end

        function computeMatFEM_FV(obj)
            % get domain ID with active Poro and Flow field
            subCells = obj.dofm.getFieldCells(obj.fields);
            nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
            iiVec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            jjVec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            Qvec = zeros((obj.mesh.nDim*obj.elements.nNodesElem)*nSubCellsByType,1);
            nDoF1 = obj.dofm.getNumDoF(obj.fields{1});
            nDoF2 = obj.dofm.getNumDoF(obj.fields{2});
            nCompPoro = obj.dofm.getDoFperEnt(obj.fields{1});
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
                    case 12 %Hexahedrons, Gauss integration
                        [N,dJWeighed] = getDerBasisFAndDet(obj.elements.hexa,el,1);
                        nG = obj.GaussPts.nNode;
                        %iN = zeros(6,8,nG); %matrix product i*N
                        B = zeros(6,8*obj.mesh.nDim,nG);
                        B(obj.elements.indB(:,2)) = N(obj.elements.indB(:,1));
                        kron = repmat([1;1;1;0;0;0],1,1,8);
                        Qs = biot*pagemtimes(B,'ctranspose',kron,'none');
                        Qs = Qs.*reshape(dJWeighed,1,1,[]);
                        Qloc = sum(Qs,3);
                        clear Qs;
                        s1 = obj.elements.nNodesElem(2)*obj.mesh.nDim;
                end
                %
                %assembly Coupling Matrix
                dofrow = getLocalDoF(obj.dofm,obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)),obj.fields{1});
                dofcol = getLocalDoF(obj.dofm,el,obj.fields{2});
                [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                iiVec(l1+1:l1+s1) = iiloc(:);
                jjVec(l1+1:l1+s1) = jjloc(:);
                Qvec(l1+1:l1+s1) = Qloc(:);
                l1 = l1+s1;
            end
            obj.Q = sparse(iiVec, jjVec, Qvec, nDoF1, nDoF2);
        end

        function stateTmp = computeRhs(obj,stateTmp,statek,dt)
            % compute Biot rhs contribute
            theta = obj.simParams.theta;
            % select active coefficients of solution vectors
            entsPoro = obj.dofm.getActiveEnts(obj.fields{1});
            entsFlow = obj.dofm.getActiveEnts(obj.fields{2});
            obj.rhs{1} = -theta*obj.Q*stateTmp.pressure(entsFlow) - ...
                (1-theta)*(obj.Q*statek.pressure(entsFlow));
            obj.rhs{2} = (obj.Q/dt)'*(stateTmp.dispCurr(entsPoro) - stateTmp.dispConv(entsPoro));
        end

        function applyDirBC(obj,field,ents,varargin)
           if strcmp(field,obj.fields{2}) && isFVTPFABased(obj.model,'Flow')
              % FV Dirichlet BC does not affect coupling blocks
              return
           else
              % call base implementation of dirichlet imposition
              applyDirBC@CouplingPhysics(obj,field,ents);
           end
        end


        function out = isLinear(obj)
            out = true;
        end


    end

end


