classdef SinglePhaseFlowFVTPFA < SinglePhaseFlow
  %SINGLEPHASEFLOW

  properties (Access = protected)
    trans
    isIntFaces
  end

  methods (Access = public)
    
    function obj = SinglePhaseFlowFVTPFA(domain)

      obj@SinglePhaseFlow(domain);

    end

    function registerSolver(obj,varargin)


      registerSolver@SinglePhaseFlow(obj,entityField.cell,varargin{:});

      %linkBoundSurf2TPFAFace(obj);

      obj.computeTrans();
      %get cells with active flow model
      flowCells = obj.domain.dofm.getActiveEntities(obj.fieldId);
      % Find internal faces (i.e. shared by two active flow cells)
      obj.isIntFaces = all(ismember(obj.grid.faces.neighbors, flowCells), 2);

      computeRhsGravTerm(obj);

    end

    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      fluid = obj.domain.materials.getFluid();
      gamma = fluid.getSpecificWeight();
      if gamma>0
        zbc = obj.grid.cells.center(:,3);
        states.potential = p + gamma*zbc;
        states.head = zbc+p/gamma;
      end
      mob = (1/fluid.getDynViscosity());
      states.flux = computeFlux(obj,p,mob,t);
      states.perm = printPermeab(obj);
      states.pressure = p;
      % states.mass = checkMassCons(obj,mob,potential);
    end

    function computeMat(obj,dt)
      % recompute elementary matrices only if the model is non-linear
      if ~isLinear(obj) || isempty(getJacobian(obj))
        mu = obj.domain.materials.getFluid().getDynViscosity();
        computeStiffMat(obj,1/mu);
        computeCapMat(obj);
      end
    end

    function computeStiffMat(obj,lw)
      % lw (the cell mobility)

      dofm = obj.domain.dofm;
      subCells = dofm.getFieldCells(obj.fieldId);
      nSubCells = length(subCells);
      % get transmissibility of internal faces
      neigh = obj.grid.faces.neighbors(obj.isIntFaces,:);
      T = lw.*obj.trans(obj.isIntFaces);
   
      neigh1 = dofm.getLocalDoF(obj.fieldId,neigh(:,1));
      neigh2 = dofm.getLocalDoF(obj.fieldId,neigh(:,2));
      sumDiagTrans = accumarray( [neigh1;neigh2], repmat(T,[2,1]), ...
        [nSubCells,1]);

      % Assemble H matrix
      nDoF = dofm.getNumbDoF(obj.fieldId);
      obj.H = sparse([neigh1; neigh2; (1:nSubCells)'],...
        [neigh2; neigh1; (1:nSubCells)'],...
        [-T; -T; sumDiagTrans], nDoF, nDoF);
    end

    function computeCapMat(obj,varargin)

      mat = obj.domain.materials;
      dofm = obj.domain.dofm;
      cellIds = dofm.getFieldCells(obj.fieldId);
      cells = obj.grid.cells;

      poroMat = zeros(cells.nTag,1);
      alphaMat = zeros(cells.nTag,1);
      beta = mat.getFluid().getFluidCompressibility();

      for m = 1:cells.nTag
        if ~ismember(m,dofm.getTargetRegions(obj.getField()))
          continue
        end

        % get regions where pressure is coupled with displacements
        coupledTags = dofm.getTargetRegions([obj.getField(),"displacements"]);
        alphaMat(m) = getRockCompressibility(obj,m,coupledTags);
        poroMat(m) = mat.getMaterial(m).PorousRock.getPorosity();

      end
      

      ctags = cells.tag(cellIds);

      % (alpha+poro*beta)
      PVal = alphaMat(ctags) + beta*poroMat(ctags);
      if ~isempty(varargin)
        % variably saturated flow model
        PVal = PVal.*varargin{1} + poroMat(ctags).*varargin{2};
      end
      PVal = PVal.*cells.volume(cellIds);
      nDoF = dofm.getNumbDoF(obj.fieldId);
      [~,~,dof] = unique(cellIds);
      obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
    end

 

    function computeRhsGravTerm(obj)
      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity'
      gamma = obj.domain.materials.getFluid().getSpecificWeight();

      if gamma > 0
        neigh = faces.neighbors(obj.isIntFaces,:);
        zVec = obj.grid.cells.center(:,3);
        zNeigh = reshape(zVec(neigh),[],2);
        obj.rhsGrav = gamma*obj.trans(obj.isIntFaces).*(zNeigh(:,1) - zNeigh(:,2));
      end
    end

    function gTerm = getRhsGravity(obj)

      fluid = obj.domain.materials.getFluid();
      lw = 1/fluid.getDynViscosity();
      
      dofm = obj.domain.dofm;
      nCells = dofm.getNumbDoF(obj.fieldId);
      neigh = obj.grid.faces.neighbors(obj.isIntFaces,:);

      gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
        -lw.*obj.rhsGrav],[nCells,1]);

      gTerm = gTerm(dofm.getActiveEntities(obj.fieldId));
      
    end

    function [ents,vals] = getBC(obj,bcId,t)
      % getBC - function to find the value and the location for the
      % boundary condition.
      %
      % Observation.:
      %  - The seepage boundary condition apply a hydrostatic pressure
      % in the boundary, and it's assume as a datum the most elavated
      % point in the domain. (For future, have a way to pass this
      % information).

      bc = obj.domain.bcs;
      mat = obj.domain.materials;
      faces = obj.grid.faces;
      surfaces = obj.grid.surfaces;
      cells = obj.grid.cells;

      bcFld = bc.getField(bcId);
      type = bc.getType(bcId);
      %
      % ents = bc.getTargetEntities(bcId);
      % vals = bc.getVals(bcId,t);

      if strcmp(type,'seepage')
        assert(isEssential(bc,bcId),"Boundary condition of type 'seepage' must be essential")
      end

      if bcFld == entityField.surface

        surfId = bc.getSourceEntities(bcId);
        faceId = surfaces.faceId(surfId);
        srcVal = bc.getSourceVals(bcId,t);
        ents = sum(faces.neighbors(faceId,:),2);
        p = getState(obj,obj.getField());

        %
        switch type
          case 'dirichlet'
            gamma = mat.getFluid().getSpecificWeight();
            mu = mat.getFluid().getDynViscosity();
            tr = obj.trans(faceId);

            dirJ = 1/mu*tr;

            dz = cells.center(ents,3) - faces.center(faceId,3);
            potential = (p(ents) - srcVal) + gamma*dz;
            q = dirJ.*potential;
            vals = [dirJ,q];
          case 'seepage'
            gamma = mat.getFluid().getSpecificWeight();
            assert(gamma>0.,'To impose Seepage boundary condition is necessary the fluid specify weight be bigger than zero!');

            zbc = faces.center(faceId,3);
            href = srcVal(1);
            v = gamma*(href-zbc);

            v(v<=0)=0.;
            mu = mat.getFluid().getDynViscosity();
            tr = obj.trans(faceId);
            dz = cells.center(ents,3) - faces.center(faceId,3);
            q = 1/mu*tr.*(p(ents) - v) + gamma*dz;
            vals = [1/mu*tr,q];
        end

      else

        ents = bc.getTargetEntities(bcId);
        vals = bc.getVals(bcId,t);

      end

    end

    function applyDirVal(obj,bcId,t)

      bcVar = obj.domain.bcs.getVariable(bcId);
      if ~strcmp(bcVar,obj.getField()) 
        return 
      end
      [bcEnts,bcVals] = getBC(obj,bcId,t);

      if size(bcVals,2)==2
        % skip BC assigned to external surfaces
        return
      end

      applyDirVal@PhysicsSolver(obj,bcId,bcEnts,bcVals);

    end

    function applyBC(obj,bcId,t)

      if ~BCapplies(obj,bcId)
        return
      end

      % helper method for FV bcs
      [bcEnts,bcVals] = getBC(obj,bcId,t);

      bcDofs = obj.domain.dofm.getLocalDoF(obj.fieldId,bcEnts);

      applyBC@PhysicsSolver(obj,bcId,bcDofs,bcVals);

    end

    function applyDirBC(obj,bcId,bcDofs,bcVals)
      % apply dirichlet BCs
      % overrides the base method implemented in PhysicsSolver
      % ents: id of constrained faces without any dof mapping applied
      % vals(:,1): Jacobian BC contrib vals(:,2): rhs BC contrib


      id = bcDofs == 0;

      bcDofs = bcDofs(~id,:);
      bcVals = bcVals(~id,:);

      % BCs imposition for finite volumes - boundary flux
      if size(bcVals,2) == 2
        assert(size(bcVals,2)==2,'Invalid matrix size for BC values');
        nDoF = obj.domain.dofm.getNumbDoF(obj.fieldId);
        bcDofsJ = nDoF*(bcDofs-1) + bcDofs;
        J = getJacobian(obj);
        obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) = ...
          J(bcDofsJ) + bcVals(:,1);
        obj.domain.rhs{obj.fieldId}(bcDofs) = ...
          obj.domain.rhs{obj.fieldId}(bcDofs) + bcVals(:,2);
      else
        applyDirBC@PhysicsSolver(obj,bcId,bcDofs);
      end
    end

    function computeTrans(obj)   % Inspired by MRST
      % Compute first the vector connecting each cell centroid to the
      % half-face

      cells = obj.grid.cells;
      faces = obj.grid.faces;

      [faceList,ptr] = getData(cells.cells2faces);
      nFPC = diff(ptr);                       % number of faces per cell

      hf2Cell = repelem((1:cells.num)',nFPC,1);
      L = faces.center(faceList,:) - cells.center(hf2Cell,:);

      % signed normal for half face transmissibility (points outside the
      % cell)
      % if a cell do not match with the first column entry of the
      % neighbors, it measn that the normal of that face points inside the
      % cell (so it has to be flipped)
      sgn = 2*(hf2Cell == faces.neighbors(faceList,1)) - 1;
      N = sgn .* faces.normal(faceList,:);

      mat = obj.domain.materials;
      KMat = zeros(cells.nTag,9);
      for i=1:cells.nTag
        KMat(i,:) = mat.getMaterial(i).PorousRock.getPermVector();
      end

      % index for vectorized L*K*N product
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];

      % half-face transmissibilities
      hT = zeros(length(hf2Cell),1);
      for k=1:length(r)
        hT = hT + L(:,r(k)) .* KMat(cells.tag(hf2Cell),k) .* N(:,c(k));
      end
      hT = hT./sum(L.*L,2);

      % compute face transmissibilities
      obj.trans = 1 ./ accumarray(faceList,1 ./ hT,[faces.num,1]);

    end


    function flux = computeFlux(obj,pot,mob,t)
      %COMPUTEFLUX - compute the flux at the faces, than accumulate
      %the value at the nodes 
      % boundary fluxes are retrieved from boundary condition value

      % Compute the nodal fluxes in the domain

      faces = obj.grid.faces;
      surfaces = obj.grid.surfaces;
      fluxFaces = zeros(faces.num,1);


      %%% process internal faces

      neigh = faces.neighbors(obj.isIntFaces,:);
      dpot = pot(neigh(:,1))-pot(neigh(:,2));
      fluxFaces(obj.isIntFaces) = mob.* obj.trans(obj.isIntFaces) * dpot;

      %%% process boundary faces

      % get boundary flux for surface bcs
      bc = obj.domain.bcs;
      bcList = bc.getBCList();
      for bcID = bcList
        if strcmp(bc.getVariable(bcID),obj.getField())
          switch bc.getField(bcID)
            case {'node','cell'}
              continue
            case 'surface'
              [~,boundFlux] = getBC(obj,bcID,t);
              surfId = bc.getSourceEntities(bcID);
              if size(vals,2) > 1 % flux associated with dirichlet bc
                boundFlux = boundFlux(:,2);
              end
          end
          faceId = surfaces.faceId(surfId);
          fluxFaces(faceId) = fluxFaces(faceId) + boundFlux;
          % add because more then one flux bc can be applied to the same
          % boundary face
        end
      end

      %%% map fluxes from face to nodes
      % provisional: nodeArea = faceArea/numNodes 
      % TO DO: use proper nodeInfluence for internal faces

      [nList,ptr] = getData(faces.connectivity);
      nNPF = diff(ptr);               % number of nodes per face
      fluxFaces = fluxFaces .* faces.normal;  % get flux direction according to computed normal

      % accumulate flux of each face on the nodes
      flux = accumarray(nList,repelem(fluxFaces./nNPF,nNPF),obj.nNodes,1);

    end


    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      pointStr = repmat(struct('name', 1, 'data', 1), 1, 1);
      pointStr(1).name = 'flux';
      pointStr(1).data = state.flux;

      cellStr = repmat(struct('name', 1, 'data', 1), 2, 1);
      cellStr(1).name = 'pressure';
      cellStr(1).data = state.pressure;
      cellStr(2).name = 'permeability';
      cellStr(2).data = state.perm;      
      if isfield(state,"potential")
        cellStr(3).name = 'potential';
        cellStr(3).data = state.potential;
        cellStr(4).name = 'piezometric head';
        cellStr(4).data = state.head;
      end
    end

    function out = isFEM(obj)
      out = false;
    end

    function out = isTPFA(obj)
      out = true;
    end

    function str = typeDiscretization(obj)
      str = "FVTPFA";
    end

  end

  methods (Access = private)

    function [lwkpt,dlwkpt] = computeMobilityBoundary(obj)
      % COMPUTEMOBILITY compute the mobility and it's derivatives
      mu = obj.domain.materials.getFluid().getDynViscosity();
      if mu==0
        lwkpt = 1;
      else
        lwkpt = 1/mu;
      end
      dlwkpt = 0;
    end

  end

end
