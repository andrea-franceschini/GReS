classdef BiotFullyCoupled < PhysicsSolver
  % Biot model subclass
  % Coupled Poromechanics with SinglePhaseFlow

  properties
    Q             % the biot coupling matrix
    R             % the fixed stress split relaxation matrix
    KdrType        % either 1,2,3 for 1D,2D or 3D bulk modulus
        
  end

  properties (Access = protected)

    % we avoid multiple inheritance and we directly create instances to the
    % single physics models that are needed

    flowSolver
    mechSolver

    % remember: all the object in the domain input are of type handle
    fldMech
    fldFlow

    flowScheme
  end

  methods (Access = public)
    function obj = BiotFullyCoupled(domain)

      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,varargin)

      dofm = obj.domain.dofm;

      default = struct('mechanicsSolver',"Poromechanics",...
                       'flowSolver',"SinglePhaseFlowFVTPFA", ...
                       'Mechanics',missing,'Flow',missing);

      input = readInput(default,varargin{:});

      % Register mechanics
      obj.mechSolver = feval(input.mechanicsSolver,obj.domain);

      % Register flow solver
      obj.flowSolver = feval(input.flowSolver,obj.domain);

      % Register mechanics solver
      registerSolver(obj.mechSolver,input.Mechanics);
      obj.fldMech = dofm.getVariableId("displacements");

      % Register flow solver
      obj.flowScheme = obj.flowSolver.typeDiscretization();
      registerSolver(obj.flowSolver,input.Flow);
      obj.fldFlow = dofm.getVariableId("pressure");

      default = struct('BulkModulus',"3D");
      params = readInput(default,input);
      blkMod = char(params.BulkModulus);
      obj.KdrType = str2double(blkMod(1));

    end

    function assembleSystem(obj,dt,varargin)

      % get Jacobian and rhs from single physics solvers
      obj.mechSolver.assembleSystem(dt);
      obj.flowSolver.assembleSystem(dt);

      % assemble coupling blocks and rhs

      % compute obj.Q
      computeMat(obj);

      % assign coupling blocks to jacobian
      obj.domain.J{obj.fldMech,obj.fldFlow} = -obj.Q;
      obj.domain.J{obj.fldFlow,obj.fldMech} = obj.Q'/dt;

      % add rhs from coupling contribution
      [rhsMech,rhsFlow] = computeRhs(obj);
      obj.domain.rhs{obj.fldMech} = obj.domain.rhs{obj.fldMech} + rhsMech;
      obj.domain.rhs{obj.fldFlow} = obj.domain.rhs{obj.fldFlow} + rhsFlow;

    end

    function computeMat(obj)

      % call method according to the discretization technique chosen
      if isempty(obj.Q) || ~isLinear(obj)
        computeMatBiot(obj)
      end
    end

    function Kdr = getDrainedBulkModulus(obj,D)

          % result is given for each gauss point

          I = zeros(6,1);
          I(1:obj.KdrType) = 1;
          DI = pagemtimes(D,I);
          Kdr = (1/obj.KdrType^2)*pagemtimes(I',DI);

    end

   function computeRelaxationMatrix(obj)

        % R = int_el b^2/Kdr*(Np'*Np)

        dofm = obj.domain.dofm;
        coordinates = obj.grid.coordinates;
        cells = obj.grid.cells;
        materials = obj.domain.materials;
        s = getState(obj);
        sOld = getStateOld(obj);
        t = s.time;

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


    function computeMatBiot(obj)
      % compute coupling matrix only where mechanics and flow are
      % active

      dofm = obj.domain.dofm;
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;

      subCells = getCoupledCells(obj);

      switch obj.flowScheme
        case "FEM"
          nEntries = sum((obj.grid.nDim)*(cells.numVerts(subCells)).^2);
        case "FVTPFA"
          nEntries = sum((obj.grid.nDim)*(cells.numVerts(subCells)));
        otherwise
          error("The discretization for the flow need to be FEM or FVTPFA.");
      end

      nDoFMech = dofm.getNumbDoF(obj.fldMech);
      nDoFFlow = dofm.getNumbDoF(obj.fldFlow);
      assembleQ = assembler(nEntries,nDoFMech,nDoFFlow);

      gaussOrd = obj.mechSolver.getGaussOrder;


      for vtkId = cells.vtkTypes

        tmp = obj.grid.getCellsByVTKId(vtkId);
        subCellsLoc = reshape(intersect(subCells,tmp,'sorted'),1,[]);
        elem = FiniteElementType.create(vtkId,obj.grid,gaussOrd);

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(subCellsLoc);
        nG = elem.getNumbGaussPts;

        for i = 1:numel(subCellsLoc)

          el = subCellsLoc(i);
          mat = obj.domain.materials.getMaterial(cells.tag(el));
          biot = mat.PorousRock.getBiotCoefficient();
 
          nodes = topol(i,:);
          coords = coordinates(nodes,:);

          % get strain matrix
          [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
          B = getStrainMatrix(elem,gradN);
          N = getBasisFinGPoints(elem);
          dofrow = getLocalDoF(dofm,obj.fldMech,nodes);

          % kronecker delta in tensor form
          kron = [1;1;1;0;0;0];
          switch obj.flowScheme
            case "FEM"
              Np = reshape(N',1,elem.nNode,nG);
              kron = repmat(kron,1,1,nG);
              iN = pagemtimes(kron,Np);
              dofcol = getLocalDoF(dofm,obj.fldFlow,nodes);
            case "FVTPFA"
              iN = repmat(kron,1,1,nG);
              dofcol = getLocalDoF(dofm,obj.fldFlow,el);
          end

          % compute local coupling matrix
          Qs = biot*pagemtimes(B,'ctranspose',iN,'none');
          Qs = Qs.*reshape(dJWeighed,1,1,[]);
          Qloc = sum(Qs,3);
          assembleQ.localAssembly(dofrow,dofcol,Qloc);

        end

      end

      obj.Q = assembleQ.sparseAssembly();

    end

    function [rhsMech,rhsFlow] = computeRhs(obj,dt)

      % retrieve State variables
      pCurr = getState(obj,"pressure");
      uCurr = getState(obj,"displacements");
      uOld = getStateOld(obj,"displacements");
      p0 = obj.getStateInit("pressure");
      dp = pCurr - p0;

      % select active coefficients of solution vectors
      entsPoro = obj.domain.dofm.getActiveEntities(obj.fldMech,1);
      entsFlow = obj.domain.dofm.getActiveEntities(obj.fldFlow,1);

      % get coupling blocks
      Qmech = getJacobian(obj,obj.fldMech,obj.fldFlow);
      Qflow = getJacobian(obj,obj.fldFlow,obj.fldMech);

      % compute rhs
      rhsMech = Qmech * dp(entsFlow);
      rhsFlow = Qflow * (uCurr(entsPoro) - uOld(entsPoro));
    end

    function initialize(obj)
      obj.mechSolver.initialize();
      obj.flowSolver.initialize();
    end

    function cells = getCoupledCells(obj)

      cellTagFlow = obj.domain.dofm.getTargetRegions(obj.fldFlow);
      cellTagMech = obj.domain.dofm.getTargetRegions(obj.fldMech);

      % find cell tag where both flow and mechanics are active
      cellTags = intersect(cellTagMech,cellTagFlow);

      cells = getEntitiesFromTags(entityField.cell,...
        obj.grid,entityField.cell,cellTags);
    end


    function flowSolv = getFlowSolver(obj)

      flowSolv = obj.flowSolver;

    end

    function mechSolv = getMechSolver(obj)

      mechSolv = obj.mechSolver;

    end

    function timeStepSetup(obj)
      obj.mechSolver.timeStepSetup();
    end


    function goBackState(obj)

      obj.mechSolver.goBackState();
    end


    function applyBC(obj,bcId,t)
      obj.flowSolver.applyBC(bcId,t);
      obj.mechSolver.applyBC(bcId,t);
    end

    function applyDirVal(obj,bcId,t)
      obj.flowSolver.applyDirVal(bcId,t);
      obj.mechSolver.applyDirVal(bcId,t);
    end


    function updateState(obj,solution)

      obj.flowSolver.updateState(solution);
      obj.mechSolver.updateState(solution);

    end

    function [cellDataBiot,pointDataBiot] = writeVTK(obj,fac,t)

      [cellDataFlow,pointDataFlow] = obj.flowSolver.writeVTK(fac,t);
      [cellDataMech,pointDataMech] = obj.mechSolver.writeVTK(fac,t);

      cellDataBiot = OutState.mergeOutFields(cellDataMech,cellDataFlow);

      clear cellDataMech celDataFlow

      pointDataBiot = OutState.mergeOutFields(pointDataMech,pointDataFlow);

      clear pointDataMech pointDataFlow

    end

    function hasConfigurationChanged = updateConfiguration(obj)

      hasConfigurationChanged = obj.mechSolver.updateConfiguration();

    end

    function writeSolution(obj,fac,tID)

      obj.flowSolver.writeSolution(fac,tID);
      obj.mechSolver.writeSolution(fac,tID);


    end


    function out = isLinear(obj)
      out = true;
    end

    function out = getFlowScheme(obj)
      out = obj.flowSolver.typeDiscretization;
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


