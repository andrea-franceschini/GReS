classdef BiotFixedStressSplit < BiotFullyCoupled
    % Biot model subclass Coupled Poromechanics with SinglePhaseFlow This
    % solver is meant for application of the fixed stress split scheme
    % Reference: Castelletto et al. (2015).
    % Accuracy and convergence properties of the fixed‐stress iterative
    % solution of two‐way coupled poromechanics
    % solved in time with backward implicit euler

    properties
        R                   % the fixed stress split relaxation matrix
        KdrType             % either 1,2,3 for 1D,2D or 3D bulk modulus
        conv                % struct with converged variables at last fixed stress iteration
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
      function obj = BiotFixedStressSplit(domain)

        obj@BiotFullyCoupled(domain);

      end


      function assembleSystem(obj,dt,var)

        % compute coupling matrix obj.Q
        computeMat(obj);

        % get Jacobian and rhs from single physics solvers
        if strcmpi(var,"pressure")
          obj.flowSolver.assembleSystem(dt);
          computeRelaxationMatrix(obj);
          rhsFlow = computeRhsFlow(bbj);
          % add the relaxation matrix 
          obj.domain.J{obj.fldFlow,obj.fldFlow} = ...
            obj.domain.J{obj.fldFlow,obj.fldFlow} + obj.R;
          obj.domain.rhs{obj.fldFlow} = obj.domain.rhs{obj.fldFlow} + rhsFlow;
        elseif strcmpi(var,"displacements")
          obj.mechSolver.assembleSystem(dt);
          rhsMech = computeRhsMech(obj);
          obj.domain.rhs{obj.fldMech} = obj.domain.rhs{obj.fldMech} + rhsMech;
        else
          error("Unknown variable for BiotFixedStressSplit solver." + ...
            "\nAvailable variables are: '%s' and '%s'",...
            obj.mechSolver.getField(),obj.flowSolver.getField())
        end

      end

      function computeMat(obj)
        
        % call method according to the discretization technique chosen
           if isempty(obj.Q) || ~isLinear(obj)
              computeMatBiot(obj)
           end
        end


        function computeMatBiot(obj)
            % compute coupling matrix only where mechanics and flow are
            % active

            subCells = getCoupledCells(obj);

            if isFEM(obj.flowSolver)
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)).^2);
            elseif isTPFA(obj.flowSolver)
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)));
            end

            [iiVec,jjVec,Qvec] = deal(zeros(nEntries,1));
            nDoF1 = obj.dofm.getNumbDoF(obj.fldMech);
            nDoF2 = obj.dofm.getNumbDoF(obj.fldFlow);
            % consider replacing the string field with an integer

            l1 = 0;
            for el=subCells'

                biot = obj.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
                elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
                nG = elem.GaussPts.nNode;
                nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));

                % get strain matrix
                [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
                iN = zeros(6,elem.nNode,nG); %matrix product i*N
                B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
                B(elem.indB(:,2)) = N(elem.indB(:,1));
                Nref = getBasisFinGPoints(elem);
                dofrow = getLocalDoF(obj.dofm,obj.fldMech,nodes);

                % kronecker delta in tensor form
                kron = [1;1;1;0;0;0];
                if isFEM(obj.flowSolver)
                  Np = reshape(Nref',1,elem.nNode,nG);
                  kron = repmat(kron,1,1,nG);
                  iN = pagemtimes(kron,Np);
                  dofcol = getLocalDoF(obj.dofm,obj.fldFlow,nodes);
                elseif isTPFA(obj.flowSolver)
                  iN = repmat(kron,1,1,nG);
                  dofcol = getLocalDoF(obj.dofm,obj.fldFlow,el);
                end

                % compute local coupling matrix
                Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
                Qs = Qs.*reshape(dJWeighed,1,1,[]);
                Qloc = sum(Qs,3);
                clear Qs;
                s1 = numel(Qloc);

                %assemble coupling Matrix
                [jjloc,iiloc] = meshgrid(dofcol,dofrow);
                iiVec(l1+1:l1+s1) = iiloc(:);
                jjVec(l1+1:l1+s1) = jjloc(:);
                Qvec(l1+1:l1+s1) = Qloc(:);
                l1 = l1+s1;
            end

            obj.Q = sparse(iiVec, jjVec, Qvec, nDoF1, nDoF2);
        end

        function computeRelaxationMatrix(obj)

          % R = int_el b^2/Kdr*(Np'*Np)

          dofm = obj.domain.dofm;

          if ~isempty(obj.R) && isLinear(obj.mechSolver)
            return
          end

          subCells = getCoupledCells(obj);

          nDoF = obj.dofm.getNumbDoF(obj.fldFlow);

          if isTPFA(obj.flowSolver)
            Rvals = zeros(nDoF,1);
            for el=subCells'
              Kdr = getDrainedBulkModulus(obj,el);
              V = obj.mesh.cellVolume(el);
              dof = dofm.getLocalDoF(obj.fldFlow,el);
              Rvals(dof) = Kdr*V;
            end
            obj.R = spdiags(Rvals,0,nDoF,nDoF);
          elseif isFEM(obj.flowSolver)
            nEntries = sum(obj.mesh.cellNumVerts(subCells).^2);
            asbR = assembler(nEntries,nDoF,nDoF);
            for el = subCells'
              Kdr = getDrainedBulkModulus(obj,el);
              nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
              elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
              dJWeighed = getDerBasisFAndDet(elem,el,3);
              N = getBasisFinGPoints(elem);
              Rloc = Kdr*(N'*diag(dJWeighed)*N);
              dof = dofm.getLocalDoF(obj.fldFlow,nodes);
              asbR.localAssembly(dof,dof,Rloc);
            end
            obj.R = asbR.sparseAssembly();
          end

        end

        function rhsMech= computeRhsMech(obj,dt)

            pCurr = getState(obj,"pressure");

            entsFlow = obj.dofm.getActiveEntities(obj.fldFlow,1);

            rhsMech = obj.Q * pCurr(entsFlow);

        end

        function [rhsMech,rhsFlow] = computeRhsFlow(obj,dt)

          % retrieve State variables
          pCurr = getState(obj,"pressure");
          uCurr = getState(obj,"displacements");
          uOld = getStateOld(obj,"displacements");
          pConv = obj.conv.pressure;
          uConv = obj.conv.displacements;

          % select active coefficients of solution vectors
          entsPoro = obj.dofm.getActiveEntities(obj.fldMech,1);
          entsFlow = obj.dofm.getActiveEntities(obj.fldFlow,1);

          % compute rhs
          rhsMech = obj.Q * pCurr(entsFlow);
          rhsFlow = obj.Q' * (uCurr(entsPoro) - uOld(entsPoro));
          % new contribution of fixed stress split algorithm
          rhsFlow = rhsFlow + obj.Q' * uConv - obj.R * pConv;
        end

        function applyBC(obj,bcId,t)
          obj.flowSolver.applyBC(bcId,t);
          obj.mechSolver.applyBC(bcId,t);
        end

        function applyDirVal(obj,bcId,t)
          obj.flowSolver.applyDirVal(bcId,t);
          obj.mechSolver.applyDirVal(bcId,t);
        end


        function updateState(obj,solution,var)
          if strmcmp(var,obj.flowSolver.getField())
            obj.flowSolver.updateState(solution);
          elseif strmcmp(var,obj.mechSolver.getField())
            obj.mechSolver.updateState(solution);
          end
        end

        function [cellDataBiot,pointDataBiot] = writeVTK(obj,t)

          [cellDataFlow,pointDataFlow] = obj.flowSolver.writeVTK(t);
          [cellDataMech,pointDataMech] = obj.mechSolver.writeVTK(t);

          cellDataBiot = OutState.mergeOutFields(cellDataMech,cellDataFlow);

          clear cellDataMech celDataFlow

          pointDataBiot = OutState.mergeOutFields(pointDataMech,pointDataFlow);

          clear pointDataMech pointDataFlow

        end

        function writeMatFile(obj,t,tID)

          obj.flowSolver.writeMatFile(t,tID);
          obj.mechSolver.writeMatFile(t,tID);


        end

        function advanceState(obj,varargin)
          if nargin == 1
            % advance state after fixed stress algorithm converged
            advanceState@BiotFullyCoupled()
          else
            % advance individual physics when newton convergence is
            % achiedeved
            varName = varargin{1};
            % save converged variable
            obj.conv.(varName) = getState(obj,varName);
          end
        end

        function getDrainedBulkModulus(obj,cTag)
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


