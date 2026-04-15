classdef SinglePhaseFlowFEM < SinglePhaseFlow
  %SINGLEPHASEFLOW

  properties (Access=protected)
    gaussOrder      % order for gauss integration (0 means the minimum required by the fem type)
  end

  methods (Access = public)
    function obj = SinglePhaseFlowFEM(domain)
      obj@SinglePhaseFlow(domain);
    end

    function registerSolver(obj,varargin)


      registerSolver@SinglePhaseFlow(obj,entityField.node,varargin{:});

      parm = readInput(struct('gaussOrder',0),varargin{:});
      obj.gaussOrder = parm.gaussOrder;

      computeRhsGravTerm(obj);

    end

    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      if gamma>0
        zbc = obj.grid.coordinates(:,3);
        states.potential = p + gamma*zbc;
        states.head = zbc+p/gamma;
      end
      states.perm = printPermeab(obj);
      states.pressure = p;
    end

    function computeMat(obj,dt)
      % recompute elementary matrices only if the model is non-linear
      if ~isLinear(obj) || isempty(getJacobian(obj))
        computeMatFEM(obj);
      end
    end

    function computeMatFEM(obj)

      % allocating vars
      dofm = obj.domain.dofm;
      materials = obj.domain.materials;
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;
      subCells = dofm.getFieldCells(obj.fieldId);
      nEntries = sum(cells.numVerts(subCells).^2);
      Ndof = dofm.getNumbDoF(obj.fieldId);
      % assembler for H and P matrix
      asbH = assembler(nEntries,Ndof,Ndof);
      asbP = assembler(nEntries,Ndof,Ndof);

      % Get the fluid compressibility
      beta = materials.getFluid().getFluidCompressibility();

      % Get the fluid dynamic viscosity
      mu = materials.getFluid().getDynViscosity();

      % get cell tags where there is coupling with the displacements
      coupledRegions = dofm.getTargetRegions([obj.getField(),"displacements"]);

      for vtkId = cells.vtkTypes

        cellList = obj.grid.getCellsByVTKId(vtkId,subCells);
        elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

        % get node topology for given vtk type
        topol = obj.grid.getCellNodes(cellList);

        for i = 1:numel(cellList)

          el = cellList(i);
          tag = cells.tag(el);
          nodes = topol(i,:);
          coords = coordinates(nodes,:);
          mat = materials.getMaterial(tag);

          permMat = mat.PorousRock.getPermMatrix();
          poro = mat.PorousRock.getPorosity();
          alpha = getRockCompressibility(obj,tag,coupledRegions);

          [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
          N = getBasisFinGPoints(elem);
          permMat = permMat/mu;
          Hs = pagemtimes(pagemtimes(gradN,'ctranspose',permMat,'none'),gradN);
          Hs = Hs.*reshape(dJWeighed,1,1,[]);
          HLoc = sum(Hs,3);
          PLoc = (alpha+poro*beta)*(N'*diag(dJWeighed)*N);
          %Getting dof associated to Flow
          dof = dofm.getLocalDoF(obj.fieldId,nodes);
          asbH.localAssembly(dof,dof,HLoc);
          asbP.localAssembly(dof,dof,PLoc);
        end

      end

      obj.H = asbH.sparseAssembly();
      obj.P = asbP.sparseAssembly();
    end


    function computeRhsGravTerm(obj)


      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity
      dofm = obj.domain.dofm;
      mat = obj.domain.materials;
      gamma = mat.getFluid().getSpecificWeight();
      coordinates = obj.grid.coordinates;
      cells = obj.grid.cells;

      if gamma > 0

        rhsTmp = zeros(dofm.getNumbDoF(obj.fieldId),1);
        subCells = dofm.getFieldCells(obj.fieldId);

        for vtkId = cells.vtkTypes

          tmp = obj.grid.getCellsByVTKId(vtkId);
          cellList = reshape(intersect(subCells,tmp,'sorted'),1,[]);
          elem = FiniteElementType.create(vtkId,obj.grid,obj.gaussOrder);

          % get node topology for given vtk type
          topol = obj.grid.getCellNodes(cellList);

          for i = 1:numel(cellList)

            % Get the material permeability
            nodes = topol(i,:);
            coords = coordinates(nodes,:);
            el = cellList(i);

            permMat = mat.getMaterial(cells.tag(el)).PorousRock.getPermMatrix();
            % why K is not divided by the dynamic viscosity?
            [gradN,dJWeighed] = getDerBasisFAndDet(elem,coords);
            fs = pagemtimes(gradN,'ctranspose',permMat(:,3),'none');
            fs = fs.*reshape(dJWeighed,1,1,[]);
            rhsLoc = sum(fs,3)*gamma;

            rhsTmp(nodes) = rhsTmp(nodes) + rhsLoc;
          end
        end

        % remove inactive components of rhs vector
        obj.rhsGrav = rhsTmp(dofm.getActiveEntities(obj.fieldId));

      end % end if
      
      end



    function rhsGrav = getRhsGravity(obj)

      rhsGrav = obj.rhsGrav;
      
    end 

    % function [ents,vals] = getBC(obj,id,t)
    %   % getBC - function to find the value and the location for the
    %   % boundary condition.
    %   %
    %   % Observation.:
    %   %  - The seepage boundary condition apply a hydrostatic pressure
    %   % in the boundary, and it's assume as a datum the most elavated
    %   % point in the domain. (For future, have a way to pass this
    %   % information).
    %   bc = obj.domain.bcs;
    % 
    %   switch bc.getCond(id)
    %     case {'node','cell'}
    %       ents = bc.getEntities(id);
    %       vals = bc.getVals(id,t);
    %     case 'surface'
    %       v = bc.getVals(id,t);
    %       ents = bc.getLoadedEntities(id);
    %       entitiesInfl = bc.getEntitiesInfluence(id);
    %       vals = entitiesInfl*v;
    %     case 'volumeforce'
    %       v = bc.getVals(id,t);
    %       ents = bc.getLoadedEntities(id);
    %       entitiesInfl = bc.getEntitiesInfluence(id);
    %       vals = entitiesInfl*v;
    %   end
    % 
    % end

    % function applyDirVal(obj,bcId,t)
    %   bcVar = obj.domain.bcs.getVariable(bcId);
    %   if ~strcmp(bcVar,obj.getField()) 
    %     return 
    %   end
    %   [bcEnts,bcVals] = getBC(obj,bcId,t);
    %   state = getState(obj);
    %   state.data.pressure(bcEnts) = bcVals;
    % end

    % function applyBC(obj,bcId,t)
    %   if ~BCapplies(obj,bcId)
    %     return
    %   end
    %   [bcDofs,bcVals] = getBC(obj,bcId,t);
    % 
    %   % Base application of a Boundary condition
    %   bcType = obj.domain.bcs.getType(bcId);
    %   switch bcType
    %     case {'dirichlet','seepage'}
    %       applyDirBC(obj,bcId,bcDofs);
    %     case {'neumann','volumeforce'}
    %       applyNeuBC(obj,bcId,bcDofs,bcVals);
    %     otherwise
    %       error("Error in %s: Boundary condition type '%s' is not " + ...
    %         "available in applyBC()",class(obj),bcType)
    %   end
    % end    

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
