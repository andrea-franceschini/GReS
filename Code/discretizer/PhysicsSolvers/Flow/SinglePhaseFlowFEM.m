classdef SinglePhaseFlowFEM < SinglePhaseFlow
  %SINGLEPHASEFLOW

  methods (Access = public)
    function obj = SinglePhaseFlowFEM(domain)
      obj@SinglePhaseFlow(domain);
    end

    function registerSolver(obj,solverInput)
      nTags = obj.mesh.nCellTag;
      
      if ~isempty(solverInput)
        targetRegions = getXMLData(solverInput,1:nTags,"targetRegions");
      else
        targetRegions = 1:nTags;
      end

      obj.dofm.registerVariable(obj.getField(),entityField.node,1,targetRegions);
      n = getNumberOfEntities(entityField.node,obj.mesh);
      obj.fieldId = obj.dofm.getVariableId(obj.getField());

      % initialize the state object with a pressure field
      obj.getState().data.(obj.getField()) = zeros(n,1);

      computeRHSGravTerm(obj);
    end

    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      gamma = obj.materials.getFluid().getSpecificWeight();
      if gamma>0
        zbc = obj.mesh.coordinates(:,3);
        states.potential = p + gamma*zbc;
        states.head = zbc+p/gamma;
      end
      states.perm = printPermeab(obj);
      states.pressure = p;
    end

    function J = computeMat(obj,dt)
      % recompute elementary matrices only if the model is non-linear
      if ~isLinear(obj) || isempty(getJacobian(obj))
        computeMatFEM(obj);
      end

      if obj.domain.simparams.isTimeDependent
        J = obj.domain.simparams.theta*obj.H + obj.P/dt;
      else
        J = obj.H;
      end
    end

    function computeMatFEM(obj)
      % dealing with input params
      subCells = obj.dofm.getFieldCells(obj.fieldId);
      nEntries = sum(obj.mesh.cellNumVerts(subCells).^2);

      [iiVec,jjVec,HVec,PVec] = deal(zeros(nEntries,1));

      % Get the fluid compressibility
      beta = obj.materials.getFluid().getFluidCompressibility();

      % Get the fluid dynamic viscosity
      mu = obj.materials.getFluid().getDynViscosity();

      l1 = 0;
      for el = subCells'
        permMat = obj.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
        poro = obj.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPorosity();
        alpha = getRockCompressibility(obj,el);
        % Compute the element matrices based on the element type
        % (tetrahedra vs. hexahedra)
        elem = getElement(obj.elements,obj.mesh.cellVTKType(el));
        [gradN,dJWeighed] = getDerBasisFAndDet(elem,el,1);
        N = getBasisFinGPoints(elem);
        permMat = permMat/mu;
        Hs = pagemtimes(pagemtimes(gradN,'ctranspose',permMat,'none'),gradN);
        Hs = Hs.*reshape(dJWeighed,1,1,[]);
        HLoc = sum(Hs,3);
        clear Hs;
        s1 = numel(HLoc);
        % Computing the P matrix contribution
        PLoc = (alpha+poro*beta)*(N'*diag(dJWeighed)*N);
        %Getting dof associated to Flow subphysic
        nodes = (obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el)));
        dof = obj.dofm.getLocalDoF(obj.fieldId,nodes);
        [jjLoc,iiLoc] = meshgrid(dof,dof);
        iiVec(l1+1:l1+s1) = iiLoc(:);
        jjVec(l1+1:l1+s1) = jjLoc(:);
        HVec(l1+1:l1+s1) = HLoc(:);
        PVec(l1+1:l1+s1) = PLoc(:);
        l1 = l1 + s1;
      end
      % renumber indices according to active nodes
      nDoF = obj.dofm.getNumbDoF(obj.fieldId);
      % Assemble H and P matrices defined as new fields of
      obj.H = sparse(iiVec, jjVec, HVec, nDoF, nDoF);
      obj.P = sparse(iiVec, jjVec, PVec, nDoF, nDoF);
    end

    function rhs = computeRhs(obj,dt)
      % Compute the residual of the flow problem

      % get pressure state
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());

      ents = obj.dofm.getActiveEntities(obj.fieldId);

      if ~obj.domain.simparams.isTimeDependent
        rhs = obj.H*p(ents);
      else
        theta = obj.domain.simparams.theta;
        rhsStiff = theta*obj.H*p(ents) + (1-theta)*obj.H*pOld(ents);
        rhsCap = (obj.P/dt)*(p(ents) - pOld(ents));
        rhs = rhsStiff + rhsCap;
      end

      %adding gravity rhs contribute
      gamma = obj.materials.getFluid().getSpecificWeight();
      if gamma > 0
        rhs = rhs + obj.rhsGrav;
      end
    end

    % TO DO: update to new FEM logic
    function computeRHSGravTerm(obj)
      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity'
      gamma = obj.materials.getFluid().getSpecificWeight();
      if gamma > 0
        rhsTmp = zeros(obj.dofm.getNumbDoF(obj.fieldId),1);
        subCells = obj.dofm.getFieldCells(obj.fieldId);
        for el = subCells'
          % Get the material permeability
          permMat = obj.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
          %             permMat = permMat/mu;
          switch obj.mesh.cellVTKType(el)
            case 10 % Tetrahedra
              N = obj.elements.tetra.getDerBasisF(el);
              rhsLoc = (N'*permMat(:,3))*obj.mesh.cellVolume(el)*gamma;
            case 12 % Hexa
              [N,dJWeighed] = obj.elements.hexa.getDerBasisFAndDet(el,1);
              fs = pagemtimes(N,'ctranspose',permMat(:,3),'none');
              fs = fs.*reshape(dJWeighed,1,1,[]);
              rhsLoc = sum(fs,3)*gamma;
          end
          entsId = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
          rhsTmp(entsId) = rhsTmp(entsId) + rhsLoc;
        end
        obj.rhsGrav = rhsTmp(obj.dofm.getActiveEnts(obj.fieldId));
      end
      % remove inactive components of rhs vector
    end

    function [ents,vals] = getBC(obj,id,t)
      % getBC - function to find the value and the location for the
      % boundary condition.
      %
      % Observation.:
      %  - The seepage boundary condition apply a hydrostatic pressure
      % in the boundary, and it's assume as a datum the most elavated
      % point in the domain. (For future, have a way to pass this
      % information).

      switch obj.bcs.getCond(id)
        case {'NodeBC','ElementBC'}
          ents = obj.bcs.getEntities(id);
          vals = obj.bcs.getVals(id,t);
        case 'SurfBC'
          v = obj.bcs.getVals(id,t);
          ents = obj.bcs.getLoadedEntities(id);
          entitiesInfl = obj.bcs.getEntitiesInfluence(id);
          vals = entitiesInfl*v;
        case 'VolumeForce'
          v = obj.bcs.getVals(id,t);
          ents = obj.bcs.getLoadedEntities(id);
          entitiesInfl = obj.bcs.getEntitiesInfluence(id);
          vals = entitiesInfl*v;
      end

    end

    function applyDirVal(obj,bcId,t)
      bcVar = obj.bcs.getVariable(bcId);
      if ~strcmp(bcVar,obj.getField()) 
        return 
      end
      [bcEnts,bcVals] = getBC(obj,bcId,t);
      state = getState(obj);
      state.data.pressure(bcEnts) = bcVals;
    end

    function applyBC(obj,bcId,t)
      if ~BCapplies(obj,bcId)
        return
      end
      [bcDofs,bcVals] = getBC(obj,bcId,t);

      % Base application of a Boundary condition
      bcType = obj.bcs.getType(bcId);
      switch bcType
        case {'Dirichlet','Seepage'}
          applyDirBC(obj,bcId,bcDofs);
        case {'Neumann','VolumeForce'}
          applyNeuBC(obj,bcId,bcDofs,bcVals);
        otherwise
          error("Error in %s: Boundary condition type '%s' is not " + ...
            "available in applyBC()",class(obj),bcType)
      end
    end    

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      cellStr = repmat(struct('name', 1, 'data', 1), 1, 1);
      cellStr(1).name = 'permeability';
      cellStr(1).data = state.perm;

      pointStr = repmat(struct('name', 1, 'data', 1), 1, 1);
      pointStr(1).name = 'pressure';
      pointStr(1).data = state.pressure;
      if isfield(state,"potential")
        pointStr(2).name = 'potential';
        pointStr(2).data = state.potential;
        pointStr(3).name = 'piezometric head';
        pointStr(3).data = state.head;
      end
    end

    function out = isFEM(obj)
      out = true;
    end

    function out = isTPFA(obj)
      out = false;
    end

    function str = typeDiscretization(obj)
      str = "FEM";
    end
  end

  methods (Access = private)

  end

end
