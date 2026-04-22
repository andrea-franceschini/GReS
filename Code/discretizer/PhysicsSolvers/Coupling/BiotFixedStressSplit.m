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

        default = struct('BulkModulus',"3D");
        params = readInput(default,input);
        blkMod = char(params.BulkModulus);

        obj.KdrType = str2double(blkMod(1));

      end


      function assembleSystem(obj,dt,var)

        % compute coupling matrix obj.Q
        computeMat(obj);

        % get Jacobian and rhs from single physics solvers
        if strcmpi(var,"pressure")
          obj.flowSolver.assembleSystem(dt);
          computeRelaxationMatrix(obj);
          rhsFlow = computeRhsFlow(obj,dt);
          % add the relaxation matrix 
          obj.domain.J{obj.fldFlow,obj.fldFlow} = ...
            obj.domain.J{obj.fldFlow,obj.fldFlow} + obj.R/dt;
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


      function computeRelaxationMatrix(obj)

        % R = int_el b^2/Kdr*(Np'*Np)

        dofm = obj.domain.dofm;
        coordinates = obj.grid.coordinates;
        cells = obj.grid.cells;
        materials = obj.domain.materials;
        s = getState(obj);
        sOld = getStateOld(obj);
        t = obj.domain.state.t;

        if ~isempty(obj.R) && isLinear(obj.mechSolver)
          return
        end

        subCells = getCoupledCells(obj);
        nDoF = dofm.getNumbDoF(obj.fldFlow);
        nEntries = sum(cells.numVerts(subCells).^2);
        asbR = assembler(nEntries,nDoF,nDoF);
        gOrd = obj.mechSolver.getGaussOrder;

        for vtkId = cells.vtkTypes

          tmp = obj.grid.getCellsByVTKId(vtkId);
          cellList = reshape(intersect(subCells,tmp,'sorted'),1,[]);
          elem = FiniteElementType.create(vtkId,obj.grid,gOrd);
          nG = elem.getNumbGaussPts;

          % get node topology for given vtk type
          topol = obj.grid.getCellNodes(cellList);

          for i = 1:numel(cellList)

            % assembly loop for homogeneous element type

            el = cellList(i);
            tag = cells.tag(el);
            nodes = topol(i,:);
            coords = coordinates(nodes,:);
            mat = materials.getMaterial(tag);

            % get drained bulk modulus
            l = obj.mechSolver.cell2stress(el,1);
            [D, ~, ~] = obj.domain.materials.updateMaterial( ...
              tag, ...
              sOld.stress(l:l+nG-1,:), ...
              s.strain(l:l+nG-1,:), ...
              [], sOld.status(l:l+nG-1,:), el, t);

            Kdr = obj.getDrainedBulkModulus(D);

            [~,dJWeighed] = getDerBasisFAndDet(elem,coords);

            if isFEM(obj.flowSolver)
              dof = dofm.getLocalDoF(obj.fldFlow,nodes);
              N = getBasisFinGPoints(elem);
            elseif isTPFA(obj.flowSolver)
              N = ones(nG,1);
              dof = dofm.getLocalDoF(obj.fldFlow,el);
            end

            biot = mat.PorousRock.getBiotCoefficient();
            k = biot^2./reshape(Kdr,1,[],1);
            Rloc = N'*diag(k.*dJWeighed)*N;
            asbR.localAssembly(dof,dof,Rloc);
          end

        end

        obj.R = asbR.sparseAssembly();
      end


      function rhsMech = computeRhsMech(obj,dt)

        pCurr = getState(obj,"pressure");
        p0 = getStateInit(obj,"pressure");

        entsFlow = obj.domain.dofm.getActiveEntities(obj.fldFlow,1);

        % remember the minus sign (look at biot)
        rhsMech = - obj.Q * (pCurr(entsFlow)-p0(entsFlow));

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


          rhsFlow = (obj.Q'/dt) * (uConv(entsPoro) - uOld(entsPoro));
          % new contribution of fixed stress split algorithm
          rhsFlow = rhsFlow  + (obj.R/dt) * (pCurr(entsFlow) - pConv(entsFlow));
        end

    

        function updateState(obj,solution,var)
          v = getState(obj,var);
          ents = obj.domain.dofm.getActiveEntities(var,1);
          v(ents) = v(ents) + solution;
          setState(obj,v,var);
        end


        function advanceState(obj,varargin)
          if nargin == 1
            % advance state after fixed stress algorithm converged
            advanceState@BiotFullyCoupled(obj)
          else
            % advance individual physics when newton convergence is
            % achiedeved
            varName = varargin{1};
            % save converged variable
            obj.conv.(varName) = getState(obj,varName);
          end
        end

        function Kdr = getDrainedBulkModulus(obj,D)

          % result is given for each gauss point

          I = zeros(6,1);
          I(1:obj.KdrType) = 1;
          DI = pagemtimes(D,I);
          Kdr = (1/obj.KdrType^2)*pagemtimes(I',DI);

        end

    end
    
    % methods (Static)
    % 
    %     function out = getField()
    %       out = [Poromechanics.getField(), SinglePhaseFlow.getField()];
    %     end
    % 
    %     function out = isSymmetric()
    %       out = false;
    %     end
    % 
    % end

end


