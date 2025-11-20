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

        obj.fldFlow = obj.dofm.getVariableId(obj.flowSolver.getField());
        obj.fldMech = obj.dofm.getVariableId(obj.mechSolver.getField());

      end

      function assembleSystem(obj,dt)

        % get Jacobian and rhs from single physics solvers
        obj.mechSolver.assembleSystem(dt);
        obj.flowSolver.assembleSystem(dt);

        % assemble coupling blocks and rhs

        % compute obj.Q
        computeMat(obj);

        % assign coupling blocks to jacobian
        obj.J{objfldMech,obj.fldFlow} = -obj.simParams.theta*obj.Q;
        obj.J{obj.fldFlow,obj.fldMech} = obj.Q'/dt;

        % add rhs from coupling contribution
        [rhsMech,rhsFlow] = computeRhs(obj);
        obj.rhs{obj.fldMech} = obj.rhs{obj.fldMech} + rhsMech;
        obj.rhs{obj.fldMech} = obj.rhs{obj.fldMech} + rhsFlow;




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

        function [rhsMech,rhsFlow] = computeRhs(obj,stateOld,dt)

            theta = obj.simParams.theta;

            % retrieve State variables
            pCurr = getState(obj,"pressure");
            pOld = getStateOld(obj,"pressure");
            uCurr = getState(obj,"displacements");
            uOld = getStateOld(obj,"displacements");

            % select active coefficients of solution vectors
            entsPoro = obj.dofm.getActiveEntities(obj.fldMech);
            entsFlow = obj.dofm.getActiveEntities(obj.fldFlow);

            % compute rhs
            rhsMech = -theta*obj.Q*pCurr(entsFlow) - (1-theta)*(obj.Q*pOld(entsFlow));
            rhsFlow = (obj.Q/dt)'*(uCurr(entsPoro) - uOld(entsPoro));
        end

        function applyBC(obj,t,bcId,bcVariable)
          obj.flowSolver.applyBC(t,bcId,bcVariable);
          obj.mechSolver.applyBC(t,bcId,bcVariable);
        end

        function applyDirVal(obj,t,bcId,bcVariable)
          obj.flowSolver.applyDirVal(t,bcId,bcVariable);
          obj.mechSolver.applyDirVal(t,bcId,bcVariable);
        end


        function updateState(obj,solution)

          obj.flowSolver.updateState(solution);
          obj.mechSolver.updateState(solution);

        end

        function [cellDataBiot,pointDataBiot] = writeVTK(obj,t)

          [cellDataFlow,pointDataFlow] = obj.flowSolver.writeVTK(t);
          [cellDataMech,pointDataMech] = obj.mechSolver.writeVTK(t);

          cellDataBiot = OutState.mergeOutFields(cellDataMech,cellDataFlow);

          clear cellDataMech celDataFlow

          pointDataBiot = OutState.mergeOutFields(pointDataMech,pointDataFlow);

          clear pointDataMech pointDataFlow

        end

        function writeMatFile(obj)

          obj.flowSolver.writeMatFile(t);
          obj.mechSolver.writeMatFile(t);


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

        function out = isSymmetric()
          out = false;
        end
        
    end

end


