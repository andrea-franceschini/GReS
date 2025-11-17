classdef newBiot < Poromechanics & SinglePhaseFlow

    properties
        biotCoefficient % example of a parameter
        J = cell(2,2) % conviently keep a cell array to separate the jacobians blocks
        rhs = cell(2,1)

    end

    properties (Constant)
      % fields array - same size of J and rhs  
      fields = ["Poromechanics","SinglePhaseFlow"];
    end

    methods (Access = public)

        function obj = Biot(inputs)
            % assign inputs
            obj.readInputs(inputs)

            % register variable fields to the dof manager
            obj.dofm.registerVariable("displacements","nodes",3,[1,2,3]);
            obj.dofm.registerVariable("pressure","cells",3,[1,2]);


        end
        
        function computeMat(obj,varargin)
           J{1,1} = computeMat@Poromechanics(obj);
           J{2,2} = computeMat@SinglePhaseFlow(obj);

           % compute the coupling blocks

        end


        function Jmat = getJacobian(obj,varargin)
            if isempty(varargin)
                Jmat = cell2mat(J);
            elseif varargin == "Poromechanics"
                Jmat = J{1,1}
                %and so on


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


