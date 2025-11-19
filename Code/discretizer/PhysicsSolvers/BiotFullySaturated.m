classdef BiotFullySaturated < PhysicsSolver
    % Biot model subclass
    % Coupled Poromechanics with SinglePhaseFlow

    properties
        Q             % the biot coupling matrix 
    end

    properties (Access = private)

      % we avoid multiple inheritance and we directly create instances to the
      % single physics models that are needed

      flowSolver
      mechSolver

      % remember: all the object in the domain input are of type handle
      fldMech
      fldFlow
    end

    methods (Access = public)
      function obj = BiotFullySaturated(domain)

        obj@PhysicsSolver(domain);

      end

      function registerSolver(obj,input)

        % setup the solver with custom input
        obj.flowSolver = SinglePhaseFlow(obj.domain);
        registerSolver(obj.flowSolver,input.(class(obj.flowSolver)));
        obj.mechSolver = Poromechanics(obj.domain);
        registerSolver(obj.mechSolver,input.(class(obj.mechSolver)));

        obj.fldMech = obj.dofm.getVariableId(obj.flowSolver.getField());
        obj.fldFlow = obj.dofm.getVariableId(obj.mechSolver.getField());

      end

      function assembleSystem(obj,dt)

        % get Jacobian and rhs from single physics solvers
        obj.mechSolver.assembleSystem(dt);
        obj.flowSolver.assembleSystem(dt);

        % assemble coupling blocks and rhs

        % compute obj.Q
        computeMat(obj);

        % assign coupling blocks
        obj.J{objfldMech,obj.fldFlow} = -obj.simParams.theta*obj.Q;
        obj.J{obj.fldFlow,obj.fldMech} = obj.Q'/dt;



      end

      function computeMat(obj)
        
        % call method according to the discretization technique chosen
           if isempty(obj.Q) || ~isLinear(obj)
              computeMatBiot(obj,dt)
           end
        end


        function computeMatBiot(obj)
            % compute coupling matrix only where mechanics and flow are
            % active

            cellTagFlow = obj.dofm.getTargetRegions(obj.fldFlow);
            cellTagMech = obj.dofm.getTargetRegions(obj.fldMech);

            cellTags = instersect(cellTagMech,cellTagFlow);

            subCells = getEntities(entityField.cell,obj.mesh,cellTags);

            if isFEM(obj.flowSolver)
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)).^2);
            elseif isFiniteVolumesTPFA(obj.flowSolver)
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)));
            end

            [iiVec,jjVec,Qvec] = deal(zeros(nEntries,1));
            nDoF1 = obj.dofm.getNumbDoF(obj.fldMech);
            nDoF2 = obj.dofm.getNumbDoF(obj.fldFlow);
            % consider replacing the string field with an integer

            l1 = 0;
            for el=subCells'
                % Get the right material stiffness for each element
                biot = obj.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
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
            entsPoro = obj.dofm.getActiveEnts(obj.fldMech);
            entsFlow = obj.dofm.getActiveEnts(obj.fldFlow);
            obj.rhs{1} = -theta*obj.Q*obj.state.data.pressure(entsFlow) - ...
                (1-theta)*(obj.Q*stateOld.data.pressure(entsFlow));
            obj.rhs{2} = (obj.Q/dt)'*(obj.state.data.dispCurr(entsPoro) - obj.state.data.dispConv(entsPoro));
        end

        function applyBC(obj,t,bcId,bcVariable)
          if strcmp(bcVariable,obj.flowSolver.getField())
          elseif strcmp(bcVariable,obj.mech.getField())
            obj.mechSolver
          end
        end


        function updateState(obj,solution)

          obj.flowSolver.updateState(solution);
          obj.mechSolver.updateState(solution);

        end

        function [cellDataBiot,pointDataBiot] = printState(obj,t)

          [cellDataFlow,pointDataFlow] = obj.flowSolver.printState(t);
          [cellDataMech,pointDataMech] = obj.flowSolver.printState(t);

          cellDataBiot = OutState.mergeOutFields(cellDataMech,cellDataFlow);
          pointDataBiot = OutState.mergeOutFields(pointDataMech,pointDataFlow);

        end

        function advanceState(obj)
          obj.mechSolver.advanceState();
          obj.flowSolver.advanceState();
        end


        function out = isLinear(obj)
            out = true;
        end
    end
    
    methods (Static)

        function out = getField()
          out = [Poromechanics.getField(), SinglePhaseFlow.getField()];
        end
    end

end


