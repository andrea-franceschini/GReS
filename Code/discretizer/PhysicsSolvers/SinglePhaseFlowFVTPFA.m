classdef SinglePhaseFlowFVTPFA < SinglePhaseFlow
  %SINGLEPHASEFLOW

  properties
    trans
    isIntFaces
  end

  methods (Access = public)
    function obj = SinglePhaseFlowFVTPFA(domain)
      obj@SinglePhaseFlow(domain);
    end

    function registerSolver(obj,solverInput)
      nTags = obj.mesh.nCellTag;

      if ~isempty(solverInput)
        targetRegions = getXMLData(solverInput,1:nTags,"targetRegions");
      else
        targetRegions = 1:nTags;
      end
      obj.dofm.registerVariable(obj.getField(),entityField.cell,1,targetRegions);
      n = getNumberOfEntities(entityField.cell,obj.mesh);
      obj.fieldId = obj.dofm.getVariableId(obj.getField());

      % initialize the state object with a pressure field
      obj.getState().data.(obj.getField()) = zeros(n,1);

      linkBoundSurf2TPFAFace(obj);

      obj.computeTrans();
      %get cells with active flow model
      flowCells = obj.dofm.getActiveEntities(obj.fieldId);
      % Find internal faces (i.e. shared by two active flow cells)
      obj.isIntFaces = all(ismember(obj.faces.faceNeighbors, flowCells), 2);

      computeRHSGravTerm(obj);
    end

    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma>0
        zbc = obj.mesh.cellCentroid(:,3);
        states.potential = p + gamma*zbc;
        states.head = zbc+p/gamma;
      end
      mob = (1/obj.materials.getFluid().getDynViscosity());
      states.flux = computeFlux(obj,p,mob,t);
      states.perm = printPermeab(obj);
      states.pressure = p;
      % states.mass = checkMassCons(obj,mob,potential);
    end

    function J = computeMat(obj,dt)
      % recompute elementary matrices only if the model is non-linear
      if ~isLinear(obj) || isempty(getJacobian(obj))
        mu = obj.materials.getFluid().getDynViscosity();
        computeStiffMat(obj,1/mu);
        computeCapMat(obj);
      end

      if obj.domain.simparams.isTimeDependent
        J = obj.domain.simparams.theta*obj.H + obj.P/dt;
      else
        J = obj.H;
      end
    end

    function computeStiffMat(obj,lw)
      % Inspired by MRST
      % subCells =
      subCells = obj.dofm.getFieldCells(obj.fieldId);
      nSubCells = length(subCells);
      %get pairs of faces that contribute to the subdomain
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      % Transmissibility of internal faces
      tmpVec = lw.*obj.trans(obj.isIntFaces);
      nneigh = length(tmpVec);
      % [~,~,reorder] = unique([neigh(:,1); neigh(:,2); subCells]);
      [~,~,reorder] = unique([neigh(:,1); neigh(:,2)]);
      neigh1 = reorder(1:nneigh);
      neigh2 = reorder(nneigh+1:2*nneigh);
      sumDiagTrans = accumarray( [neigh1;neigh2], repmat(tmpVec,[2,1]), ...
        [nSubCells,1]);
      % Assemble H matrix
      nDoF = obj.dofm.getNumbDoF(obj.fieldId);
      obj.H = sparse([neigh1; neigh2; (1:nSubCells)'],...
        [neigh2; neigh1; (1:nSubCells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], nDoF, nDoF);
    end

    function computeCapMat(obj,varargin)
      subCells = obj.dofm.getFieldCells(obj.fieldId);
      nSubCells = length(subCells);
      poroMat = zeros(nSubCells,1);
      alphaMat = zeros(nSubCells,1);
      beta = obj.materials.getFluid().getFluidCompressibility();
      for m = 1:obj.mesh.nCellTag
        if ~ismember(m,obj.dofm.getTargetRegions([obj.getField(),"displacements"]))
          % compute alpha only if there's no coupling in the
          % subdomain
          alphaMat(m) = obj.materials.getMaterial(m).ConstLaw.getRockCompressibility();
        end
        poroMat(m) = obj.materials.getMaterial(m).PorousRock.getPorosity();
      end
      % (alpha+poro*beta)
      PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));
      if ~isempty(varargin)
        % variably saturated flow model
        PVal = PVal.*varargin{1} + poroMat(obj.mesh.cellTag(subCells)).*varargin{2};
      end
      PVal = PVal.*obj.mesh.cellVolume(subCells);
      nDoF = obj.dofm.getNumbDoF(obj.fieldId);
      [~,~,dof] = unique(subCells);
      obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
    end

    function rhs = computeRhs(obj,dt)
      % Compute the residual of the flow problem

      % get pressure state
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());

      lw = 1/obj.materials.getFluid().getDynViscosity();
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
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma > 0
        rhs = rhs + finalizeRHSGravTerm(obj,lw);
      end

    end

    function computeRHSGravTerm(obj)
      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity'
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma > 0
        neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
        zVec = obj.mesh.cellCentroid(:,3);
        zNeigh = zVec(neigh);
        obj.rhsGrav = gamma*obj.trans(obj.isIntFaces).*(zNeigh(:,1) - zNeigh(:,2));
      end
      % remove inactive components of rhs vector
    end

    function gTerm = finalizeRHSGravTerm(obj,lw)
      nCells = obj.dofm.getNumbDoF(obj.fieldId);
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
        -lw.*obj.rhsGrav],[nCells,1]);
      gTerm = gTerm(obj.dofm.getActiveEntities(obj.fieldId));
    end

    function [dof,vals] = getBC(obj,id,t)
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

          faceID = obj.bcs.getEntities(id);
          ents = sum(obj.faces.faceNeighbors(faceID,:),2);

          p = getState(obj,"pressure");

          % [ents,~,ind] = unique(ents);
          % % % [faceID, faceOrder] = sort(obj.bcs.getEntities(id));
          % % % ents = sum(obj.faces.faceNeighbors(faceID,:),2);
          % % % v(faceOrder,1) = obj.bcs.getVals(id,t);

          switch obj.bcs.getType(id)

            case 'Neumann'
              vals = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;

            case 'Dirichlet'
              gamma = obj.materials.getFluid().getFluidSpecWeight();
              mu = obj.materials.getFluid().getDynViscosity();
              tr = obj.trans(faceID);

              % q = 1/mu*tr.*((obj.state.data.pressure(ents) - v)...
              %    + gamma*(obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3)));
              % vals = [1/mu*tr,accumarray(ind,q)]; % {JacobianVal,rhsVal]

              dirJ = 1/mu*tr;
              % % press = obj.state.data.pressure(ents) - v;
              % % gravT =  gamma*(obj.mesh.cellCentroid(ents,3) ...
              % %   - obj.faces.faceCentroid(faceID,3));
              dz = obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3);
              potential = (p(ents) - v) + gamma*dz;
              q = dirJ.*potential;
              vals = [dirJ,q];

            case 'Seepage'
              gamma = obj.materials.getFluid().getFluidSpecWeight();
              assert(gamma>0.,'To impose Seepage boundary condition is necessary the fluid specify weight be bigger than zero!');

              zbc = obj.faces.faceCentroid(faceID,3);
              href = v(1);
              v = gamma*(href-zbc);

              v(v<=0)=0.;
              mu = obj.materials.getFluid().getDynViscosity();
              tr = obj.trans(faceID);
              dz = obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3);
              q = 1/mu*tr.*(p(ents) - v) + gamma*dz;
              vals = [1/mu*tr,q];
          end


        case 'VolumeForce'
          v = obj.bcs.getVals(id,t);
          ents = obj.bcs.getEntities(id);
          vals = v.*obj.mesh.cellVolume(ents);
      end

      dof = obj.dofm.getLocalDoF(obj.fieldId,ents);
    end

    function applyDirVal(obj,bcId,t)
      return
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
          applyDirBC(obj,bcId,bcDofs,bcVals);
        case {'Neumann','VolumeForce'}
          applyNeuBC(obj,bcId,bcDofs,bcVals);
        otherwise
          error("Error in %s: Boundary condition type '%s' is not " + ...
            "available in applyBC()",class(obj),bcType)
      end
    end

    function applyDirBC(obj,bcId,bcDofs,bcVals)
      % apply Dirichlet BCs
      % ents: id of constrained faces without any dof mapping applied
      % vals(:,1): Jacobian BC contrib vals(:,2): rhs BC contrib

      % BCs imposition for finite volumes - boundary flux
      assert(size(bcVals,2)==2,'Invalid matrix size for BC values');
      nDoF = obj.dofm.getNumbDoF(obj.fieldId);
      bcDofsJ = nDoF*(bcDofs-1) + bcDofs;
      J = getJacobian(obj);
      obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) = ...
        J(bcDofsJ) + bcVals(:,1);
      obj.domain.rhs{obj.fieldId}(bcDofs) = ...
        obj.domain.rhs{obj.fieldId}(bcDofs) + bcVals(:,2);
    end

    function computeTrans(obj)   % Inspired by MRST
      % Compute first the vector connecting each cell centroid to the
      % half-face
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];
      hf2Cell = repelem((1:obj.mesh.nCells)',diff(obj.faces.mapF2E));
      L = obj.faces.faceCentroid(obj.faces.faces2Elements(:,1),:) - obj.mesh.cellCentroid(hf2Cell,:);
      sgn = 2*(hf2Cell == obj.faces.faceNeighbors(obj.faces.faces2Elements(:,1))) - 1;
      N = bsxfun(@times,sgn,obj.faces.faceNormal(obj.faces.faces2Elements(:,1),:));
      KMat = zeros(obj.mesh.nCellTag,9);
      for i=1:obj.mesh.nCellTag
        KMat(i,:) = obj.materials.getMaterial(i).PorousRock.getPermVector();
      end
      hT = zeros(length(hf2Cell),1);
      for k=1:length(r)
        hT = hT + L(:,r(k)) .* KMat(obj.mesh.cellTag(hf2Cell),k) .* N(:,c(k));
      end
      hT = hT./sum(L.*L,2);
      %       mu = obj.materials.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
      %       hT = hT/mu;
      obj.trans = 1 ./ accumarray(obj.faces.faces2Elements(:,1),1 ./ hT,[obj.faces.nFaces,1]);
    end

    function trans = computeTransCell(obj,ncell)
      % Compute first the vector connecting each cell centroid to the
      % half-face
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];
      nrep = diff(obj.faces.mapF2E);
      hf2Cell = repelem(ncell.gtCell+1,nrep(ncell.gtCell+1));
      L = obj.faces.faceCentroid(ncell.face,:) - obj.elements.cellCentroid(hf2Cell,:);
      % sgn = 2*(hf2Cell == obj.faces.faceNeighbors(ncell.face)) - 1;
      sgn(1:6)=-1;
      N = bsxfun(@times,sgn,obj.faces.faceNormal(ncell.face,:)')';
      KMat = zeros(obj.mesh.nCellTag,9);
      for i=1:obj.mesh.nCellTag
        KMat(i,:) = obj.materials.getMaterial(i).PorousRock.getPermVector();
      end
      hT = zeros(length(hf2Cell),1);
      for k=1:length(r)
        hT = hT + L(:,r(k)) .* KMat(obj.mesh.cellTag(hf2Cell),k) .* N(:,c(k));
      end
      trans = hT./sum(L.*L,2);
    end

    function flux = computeFlux(obj,pot,mob,t)
      %COMPUTEFLUX - compute the flux at the faces, than accumulate
      %the value at the nodes (The contribution of the boundary is done
      % in another function).

      % Compute the fluxes inside the domain.
      nnodesBfaces = diff(obj.faces.mapN2F);
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      % isIntNode = repelem(obj.isIntFaces,nnodesBfaces);

      fluxFaces = zeros(obj.faces.nFaces,1);
      fluxFaces(obj.isIntFaces) = pot(neigh(:,1))-pot(neigh(:,2));
      fluxFaces(obj.isIntFaces) = mob.*obj.trans(obj.isIntFaces).*fluxFaces(obj.isIntFaces);
      % fluxFaces = pot(neigh(:,1))-pot(neigh(:,2));
      % fluxFaces = mob.*obj.trans(obj.isIntFaces).*fluxFaces;

      areaSq = vecnorm(obj.faces.faceNormal,2,2);
      faceUnit = obj.faces.faceNormal./areaSq;
      areaSq = areaSq.*nnodesBfaces;
      fluxFaces = fluxFaces./areaSq;
      fluxFaces = fluxFaces.*faceUnit;
      % flux at the faces
      % fluxFaces = fluxFaces./areaSq(obj.isIntFaces);
      % fluxFaces = fluxFaces.*faceUnit(obj.isIntFaces,:);  % flux at the faces

      % Contribution
      fluxFaces = repelem(fluxFaces,nnodesBfaces,1);
      axis = ones(length(obj.faces.nodes2Faces),1);
      flux = accumarray([[obj.faces.nodes2Faces axis]; ...
        [obj.faces.nodes2Faces 2*axis]; ...
        [obj.faces.nodes2Faces 3*axis]], fluxFaces(:));

      % fluxFaces = repelem(fluxFaces,nnodesBfaces(obj.isIntFaces),1);
      % axis = ones(size(fluxFaces,1),1);
      % flux = accumarray([[obj.faces.nodes2Faces(isIntNode) axis]; ...
      %    [obj.faces.nodes2Faces(isIntNode) 2*axis]; ...
      %    [obj.faces.nodes2Faces(isIntNode) 3*axis]], fluxFaces(:));

      % Compute the fluxes at the boundary of the domain.
      Node2Face = repelem((1:obj.faces.nFaces)',nnodesBfaces);
      sgn = 2*(obj.faces.faceNeighbors(:,1)==0) - 1;

      areaSq = vecnorm(obj.faces.faceNormal,2,2);
      faceUnit = obj.faces.faceNormal./areaSq;
      areaSq = areaSq.*nnodesBfaces;

      % add boundary condition
      bcList = obj.bcs.db.keys;
      for bc = string(bcList)
        if strcmp(obj.bcs.getVariable(bc),obj.getField())
          [~,vals] = getBC(obj,bc,t);
          [ents, ~] = sort(obj.bcs.getEntities(bc));
          switch obj.bcs.getCond(bc)
            case {'NodeBC','ElementBC'}
            case 'SurfBC'
              if strcmp(obj.bcs.getType(bc),'Dirichlet') || strcmp(obj.bcs.getType(bc),'Seepage')
                vals=vals(:,2);
              end
              dir = sgn(ents).*faceUnit(ents,:);
              vals = vals(:)./areaSq(ents).*dir;
              vals = repelem(vals,nnodesBfaces(ents),1);
            case 'VolumeForce'
              facesBcell = diff(obj.faces.mapF2E);

              % Find the faces to distribute the contribution.
              vals = vals(:)./facesBcell(ents);
              vals = repelem(vals,facesBcell(ents),1);

              hf2Cell = repelem((1:obj.mesh.nCells)',facesBcell);
              ents = obj.faces.faces2Elements(hf2Cell == ents,1);

              vals = sgn(ents).*vals./areaSq(ents).*faceUnit(ents,:);
              vals = -repelem(vals,nnodesBfaces(ents),1);
          end

          nodes = obj.faces.nodes2Faces(ismember(Node2Face,ents));
          [loc,~,pos] = unique(nodes);
          axis = ones(length(nodes),1);
          fluxB = accumarray([[pos axis]; [pos 2*axis]; ...
            [pos 3*axis]], vals(:));
          flux(loc,:)=fluxB;
        end
      end
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

    % function mass = checkMassCons(obj,mob,pot)
    %   %CHECKMASSCONS - check the mass conservation in all elements.
    %   mass = zeros(obj.mesh.nCells,1);
    %   if isTPFA(obj)
    %     neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
    %     sgn = 2*((obj.faces.faces2Elements(:,2)==1) +(obj.faces.faces2Elements(:,2)==3)+(obj.faces.faces2Elements(:,2)==5)) - 1;
    %
    %     fluxFaces = zeros(obj.faces.nFaces,1);
    %     fluxFaces(obj.isIntFaces) = pot(neigh(:,1))-pot(neigh(:,2));
    %     fluxFaces(obj.isIntFaces) = mob.*obj.trans(obj.isIntFaces).*fluxFaces(obj.isIntFaces);
    %
    %     % Contribution
    %     massFace = sgn.*fluxFaces(obj.faces.faces2Elements(:,1));
    %     elm = repelem(1:obj.mesh.nCells,diff(obj.faces.mapF2E));
    %     mass = accumarray(elm',massFace);
    %   end
    % end
  end

  methods (Access = private)

    function [lwkpt,dlwkpt] = computeMobilityBoundary(obj)
      % COMPUTEMOBILITY compute the mobility and it's derivatives
      mu = obj.materials.getFluid().getDynViscosity();
      if mu==0
        lwkpt = 1;
      else
        lwkpt = 1/mu;
      end
      dlwkpt = 0;
    end

  end

end
