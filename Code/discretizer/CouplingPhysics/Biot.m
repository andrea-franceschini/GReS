classdef Biot < CouplingPhysics
    % Biot model subclass
    % Coupled poromechanics with Flow models
    % field1: Poromechanics
    % field2: SPFlow

    properties
        Q
        flowModel   % Model used for Flow
        flowScheme  % Discretization scheme used for flow
    end

    properties (Constant)
      fields = ["Poromechanics","SinglePhaseFlow"];
    end

    methods (Access = public)
        function obj = Biot(symmod,params,dofManager,grid,mat,bc,state)
            obj@CouplingPhysics(symmod,params,dofManager,grid,mat,bc,state);
            if isSinglePhaseFlow(obj.model)
                obj.flowModel = 'SinglePhaseFlow';
            elseif isVariabSatFlow(obj.model)
                obj.flowModel = 'VaraiblySaturatedFlow';
            end
            obj.flowScheme = obj.dofm.getScheme(obj.fields(2));
            %
        end

        function computeMat(obj,varargin)
           dt = varargin{2};
           % call method according to the discretization technique chosen
           if isempty(obj.J{1}) || ~isLinear(obj)
              computeMatBiot(obj)
           end
           % J{1}: momentum balance   J{2}: mass balance equations
           obj.J{1} = -obj.simParams.theta*obj.Q;
           obj.J{2} = obj.Q'/dt;
        end


        function computeMatBiot(obj)
            % get domain ID with active Poro and Flow field
            subCells = obj.dofm.getFieldCells(obj.fields);
            switch obj.flowScheme
              case 'FEM'
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)).^2);
              case 'FV'
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)));
            end
            [iiVec,jjVec,Qvec] = deal(zeros(nEntries,1));
            nDoF1 = obj.dofm.getNumDoF(obj.fields(1));
            nDoF2 = obj.dofm.getNumDoF(obj.fields(2));
            %
            l1 = 0;
            for el=subCells'
                % Get the right material stiffness for each element
                biot = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
                elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
                nG = elem.GaussPts.nNode;
                nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
                [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
                iN = zeros(6,elem.nNode,nG); %matrix product i*N
                B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
                B(elem.indB(:,2)) = N(elem.indB(:,1));
                Nref = getBasisFinGPoints(elem);
                % kronecker delta in tensor form
                dofrow = getLocalDoF(obj.dofm,nodes,obj.fldId(1));
                kron = [1;1;1;0;0;0];
                switch obj.flowScheme
                  case 'FEM'
                    Np = reshape(Nref',1,elem.nNode,nG);
                    kron = repmat(kron,1,1,nG);
                    iN = pagemtimes(kron,Np);
                    dofcol = getLocalDoF(obj.dofm,nodes,obj.fldId(2));
                  case 'FV'
                    iN = repmat(kron,1,1,nG);
                    dofcol = getLocalDoF(obj.dofm,el,obj.fldId(2));
                end
                Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
                Qs = Qs.*reshape(dJWeighed,1,1,[]);
                Qloc = sum(Qs,3);
                clear Qs;
                s1 = numel(Qloc);
                %assembly Coupling Matrix
                [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                iiVec(l1+1:l1+s1) = iiloc(:);
                jjVec(l1+1:l1+s1) = jjloc(:);
                Qvec(l1+1:l1+s1) = Qloc(:);
                l1 = l1+s1;
            end
            obj.Q = sparse(iiVec, jjVec, Qvec, nDoF1, nDoF2);
        end

        function computeRhs(obj,stateOld,dt)
            % compute Biot rhs contribute
            theta = obj.simParams.theta;
            % select active coefficients of solution vectors
            entsPoro = obj.dofm.getActiveEnts(obj.fields(1));
            entsFlow = obj.dofm.getActiveEnts(obj.fields(2));
            obj.rhs{1} = -theta*obj.Q*obj.state.data.pressure(entsFlow) - ...
                (1-theta)*(obj.Q*stateOld.data.pressure(entsFlow));
            obj.rhs{2} = (obj.Q/dt)'*(obj.state.data.dispCurr(entsPoro) - obj.state.data.dispConv(entsPoro));
        end

        function applyDirBC(obj,field,ents,varargin)
           if strcmp(field,obj.fields(2)) && isFVTPFABased(obj.model,'Flow')
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
    
    methods (Static)

        function out = getField()
          out = Biot.fields;
        end
    end

end


