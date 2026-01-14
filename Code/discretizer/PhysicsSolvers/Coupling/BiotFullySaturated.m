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

    flowScheme
  end

  methods (Access = public)
    function obj = BiotFullySaturated(domain)

      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,input)

      % Register mechanics
      obj.mechSolver = Poromechanics(obj.domain);
      registerSolver(obj.mechSolver,input.(class(obj.mechSolver)));
      obj.fldMech = obj.dofm.getVariableId(obj.mechSolver.getField());

      % setup the solver with custom input
      if isfield(input,"SinglePhaseFlowFEM")
        obj.flowSolver = SinglePhaseFlowFEM(obj.domain);
      end
      if isfield(input,"SinglePhaseFlowFVTPFA")
        obj.flowSolver = SinglePhaseFlowFVTPFA(obj.domain);
      end

      % Register fluids
      obj.flowScheme = obj.flowSolver.typeDiscretization();
      registerSolver(obj.flowSolver,input.(class(obj.flowSolver)));
      obj.fldFlow = obj.dofm.getVariableId(obj.flowSolver.getField());

    end

    function assembleSystem(obj,dt)

      % get Jacobian and rhs from single physics solvers
      obj.mechSolver.assembleSystem(dt);
      obj.flowSolver.assembleSystem(dt);

      % assemble coupling blocks and rhs

      % compute obj.Q
      computeMat(obj);

      % assign coupling blocks to jacobian
      obj.domain.J{obj.fldMech,obj.fldFlow} = -obj.domain.simparams.theta*obj.Q;
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


    function computeMatBiot(obj)
      % compute coupling matrix only where mechanics and flow are
      % active

      cellTagFlow = obj.dofm.getTargetRegions(obj.fldFlow);
      cellTagMech = obj.dofm.getTargetRegions(obj.fldMech);

      % find cell tag where both flow and mechanics are active
      cellTags = intersect(cellTagMech,cellTagFlow);

      subCells = getEntities(entityField.cell,obj.mesh,cellTags);

      switch obj.flowScheme
        case "FEM"
          nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)).^2);
        case "FVTPFA"
          nEntries = sum((obj.mesh.nDim)*(obj.mesh.cellNumVerts(subCells)));
        otherwise
          error("The discretization for the flow need to be FEM or FVTPFA.");
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
        switch obj.flowScheme
          case "FEM"
            Np = reshape(Nref',1,elem.nNode,nG);
            kron = repmat(kron,1,1,nG);
            iN = pagemtimes(kron,Np);
            dofcol = getLocalDoF(obj.dofm,obj.fldFlow,nodes);
          case "FVTPFA"
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

    function [rhsMech,rhsFlow] = computeRhs(obj,dt)

      % retrieve State variables
      pCurr = getState(obj,"pressure");
      pOld = getStateOld(obj,"pressure");
      uCurr = getState(obj,"displacements");
      uOld = getStateOld(obj,"displacements");

      % select active coefficients of solution vectors
      entsPoro = obj.dofm.getActiveEntities(obj.fldMech,1);
      entsFlow = obj.dofm.getActiveEntities(obj.fldFlow,1);

      % get coupling blocks
      Qmech = getJacobian(obj,obj.fldMech,obj.fldFlow);
      Qflow = getJacobian(obj,obj.fldFlow,obj.fldMech);

      % compute rhs
      theta = obj.domain.simparams.theta;
      rhsMech = Qmech * (pCurr(entsFlow) + (1/theta-1)*pOld(entsFlow));
      rhsFlow = Qflow * (uCurr(entsPoro) - uOld(entsPoro));
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
          [cellDataMech,pointDataMech] = obj.mechSolver.writeVTK(fac);

      cellDataBiot = OutState.mergeOutFields(cellDataMech,cellDataFlow);

      clear cellDataMech celDataFlow

      pointDataBiot = OutState.mergeOutFields(pointDataMech,pointDataFlow);

      clear pointDataMech pointDataFlow

    end

    function writeMatFile(obj,fac,tID)

      obj.flowSolver.writeMatFile(fac,tID);
      obj.mechSolver.writeMatFile(fac,tID);


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


