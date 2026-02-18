classdef BiotFixedStressSplit < BiotFullyCoupled
    % Biot model subclass of Coupled Poromechanics with SinglePhaseFlow 
    % This solver is meant for application of the fixed stress split scheme
    % Reference: Castelletto et al. (2015). Accuracy and convergence
    % properties of the fixed‐stress iterative solution of two‐way coupled
    % poromechanics solved in time with backward implicit euler

    properties
        R              % the fixed stress split relaxation matrix
        KdrType        % either 1,2,3 for 1D,2D or 3D bulk modulus
        conv           % struct with converged variables at last fixed stress iteration
    end


    methods (Access = public)
      function obj = BiotFixedStressSplit(domain)

        obj@BiotFullyCoupled(domain);

      end

      function registerSolver(obj,input)

        registerSolver@BiotFullyCoupled(obj,input);

        kdr = char(getXMLData(input,"3D","BulkModulus"));
        obj.KdrType = str2double(kdr(1));

        % obj.conv.(obj.flowSolver.getField) = getState(obj,obj.flowSolver.getField);
        % obj.conv.(obj.mechSolver.getField) = getState(obj,obj.mechSolver.getField);

      end


      function assembleSystem(obj,dt,var)

        % compute coupling matrix obj.Q
        computeMat(obj);

        % get Jacobian and rhs from single physics solvers
        if strcmpi(var,"pressure")
          obj.flowSolver.assembleSystem(dt);
          computeRelaxationMatrix(obj);
          rhsFlow = computeRhsFlow(obj);
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

            dofm = obj.domain.dofm;

            subCells = getCoupledCells(obj);

            if isFEM(obj.flowSolver)
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)).^2);
            elseif isTPFA(obj.flowSolver)
                nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)));
            end

            [iiVec,jjVec,Qvec] = deal(zeros(nEntries,1));
            nDoF1 = dofm.getNumbDoF(obj.fldMech);
            nDoF2 = dofm.getNumbDoF(obj.fldFlow);
            % consider replacing the string field with an integer

            l1 = 0;
            for el=subCells'

                biot = obj.domain.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getBiotCoefficient();
                elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
                nG = elem.GaussPts.nNode;
                nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));

                % get strain matrix
                [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
                iN = zeros(6,elem.nNode,nG); %matrix product i*N
                B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
                B(elem.indB(:,2)) = N(elem.indB(:,1));
                Nref = getBasisFinGPoints(elem);
                dofrow = getLocalDoF(dofm,obj.fldMech,nodes);

                % kronecker delta in tensor form
                kron = [1;1;1;0;0;0];
                if isFEM(obj.flowSolver)
                  Np = reshape(Nref',1,elem.nNode,nG);
                  kron = repmat(kron,1,1,nG);
                  iN = pagemtimes(kron,Np);
                  dofcol = getLocalDoF(dofm,obj.fldFlow,nodes);
                elseif isTPFA(obj.flowSolver)
                  iN = repmat(kron,1,1,nG);
                  dofcol = getLocalDoF(dofm,obj.fldFlow,el);
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

        % function computeRelaxationMatrix(obj)
        % 
        %   % R = int_el b^2/Kdr*(Np'*Np)
        % 
        %   dofm = obj.domain.dofm;
        % 
        %   if ~isempty(obj.R) && isLinear(obj.mechSolver)
        %     return
        %   end
        % 
        %   subCells = getCoupledCells(obj);
        % 
        %   nDoF = dofm.getNumbDoF(obj.fldFlow);
        % 
        %   if isTPFA(obj.flowSolver)
        %     Rvals = zeros(nDoF,1);
        %     for el=subCells'
        %       Kdr = getDrainedBulkModulus(obj,el);
        %       V = obj.mesh.cellVolume(el);
        %       dof = dofm.getLocalDoF(obj.fldFlow,el);
        %       Rvals(dof) = Kdr*V;
        %     end
        %     obj.R = spdiags(Rvals,0,nDoF,nDoF);
        %   elseif isFEM(obj.flowSolver)
        %     nEntries = sum(obj.mesh.cellNumVerts(subCells).^2);
        %     asbR = assembler(nEntries,nDoF,nDoF);
        %     for el = subCells'
        %       Kdr = getDrainedBulkModulus(obj,el);
        %       nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
        %       elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
        %       dJWeighed = getDerBasisFAndDet(elem,el,3);
        %       N = getBasisFinGPoints(elem);
        %       Rloc = N'*diag(Kdr.*dJWeighed)*N;
        %       dof = dofm.getLocalDoF(obj.fldFlow,nodes);
        %       asbR.localAssembly(dof,dof,Rloc);
        %     end
        %     obj.R = asbR.sparseAssembly();
        %   end
        %
        % end


        function computeRelaxationMatrix(obj)

          % R = int_el b^2/Kdr*(Np'*Np)

          dofm = obj.domain.dofm;

          if ~isempty(obj.R) && isLinear(obj.mechSolver)
            return
          end

          subCells = getCoupledCells(obj);

          nDoF = dofm.getNumbDoF(obj.fldFlow);

          nEntries = sum(obj.mesh.cellNumVerts(subCells).^2);
          asbR = assembler(nEntries,nDoF,nDoF);


          for el = subCells'
            % get bulk modulus
            Kdr = getDrainedBulkModulus(obj,el);
            mat = obj.domain.materials.getMaterial(obj.mesh.cellTag(el));
            elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
            dJWeighed = getDerBasisFAndDet(elem,el,3);
            if isFEM(obj.flowSolver)
              nodes = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
              dof = dofm.getLocalDoF(obj.fldFlow,nodes);
              N = getBasisFinGPoints(elem);
            elseif isTPFA(obj.flowSolver)
              nG = elem.GaussPts.nNode;
              N = ones(nG,1);
              dof = dofm.getLocalDoF(obj.fldFlow,el);
            end
            biot = mat.PorousRock.getBiotCoefficient();
            k = biot^2./reshape(Kdr,1,[],1);
            Rloc = N'*diag(k.*dJWeighed)*N;
            asbR.localAssembly(dof,dof,Rloc);
          end

          obj.R = asbR.sparseAssembly();
        end


        function rhsMech= computeRhsMech(obj,dt)

          pCurr = getState(obj,"pressure");

          entsFlow = obj.domain.dofm.getActiveEntities(obj.fldFlow,1);

          rhsMech = obj.Q * pCurr(entsFlow);

        end

        function rhsFlow = computeRhsFlow(obj,dt)

          % retrieve State variables
          pCurr = getState(obj,"pressure");
          uOld = getStateOld(obj,"displacements");
          pConv = obj.conv.pressure;
          uConv = obj.conv.displacements;

          % select active coefficients of solution vectors
          entsPoro = obj.domain.dofm.getActiveEntities(obj.fldMech,1);
          entsFlow = obj.domain.dofm.getActiveEntities(obj.fldFlow,1);


          rhsFlow = obj.Q' * (uConv(entsPoro) - uOld(entsPoro));
          % new contribution of fixed stress split algorithm
          rhsFlow = rhsFlow  + obj.R * (pCurr(entsFlow) - pConv(entsFlow));
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
          if strcmp(var,obj.flowSolver.getField())
            ents = obj.domain.dofm.getActiveEntities(obj.fldFlow,1);
            state = getState(obj);
            state.data.pressure(ents) = state.data.pressure(ents) + solution;
          elseif strcmp(var,obj.mechSolver.getField())
            ents = obj.domain.dofm.getActiveEntities(obj.fldMech,1);
            state = getState(obj);
            state.data.displacements(ents) = state.data.displacements(ents) + solution;
            obj.mechSolver.updateState();
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

        function Kdr = getDrainedBulkModulus(obj,elID)

          % result is given for each gauss point

          s = getState(obj);
          sOld = getStateOld(obj);

          vtkId = obj.mesh.cellVTKType(elID);
          elem = getElement(obj.elements,vtkId);
          nG = elem.GaussPts.nNode;

          l = obj.mechSolver.cell2stress(elID);

          % get constitutive matrix
          [D, ~, ~] = obj.domain.materials.updateMaterial( ...
            obj.mesh.cellTag(elID), ...
            sOld.data.stress(l+1:l+nG,:), ...
            s.data.strain(l+1:l+nG,:), ...
            [], sOld.data.status(l+1:l+nG,:), elID, s.t);

          I = zeros(6,1);
          I(1:obj.KdrType) = 1;
          DI = pagemtimes(D,I);
          Kdr = (1/obj.KdrType^2)*pagemtimes(I',DI);

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


