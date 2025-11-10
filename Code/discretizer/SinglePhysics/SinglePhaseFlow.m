classdef SinglePhaseFlow < SinglePhysics
  %POROMECHANICS
  % Subclass of Discretizer
  % Implement Poromechanics methods to assemble the stiffness matrix and
  % the residual contribution

  properties
    trans
    isIntFaces
    rhsGrav         % gravity contribution to rhs
    H
    P
  end

  properties (Constant)
    field = 'SinglePhaseFlow'
  end

  methods (Access = public)
    function obj = SinglePhaseFlow(symmod,params,dofManager,grid,mat,bc,state)
      obj@SinglePhysics(symmod,params,dofManager,grid,mat,bc,state);
      if obj.model.isFVTPFABased('Flow')
        obj.computeTrans;
        %get cells with active flow model
        flowCells = obj.dofm.getFieldCells(obj.getField());
        % Find internal faces (i.e. shared by two cells with active
        % flow)
        obj.isIntFaces = all(ismember(obj.faces.faceNeighbors, flowCells), 2);
      end
      computeRHSGravTerm(obj);
    end

    function initState(obj)
      if isFEMBased(obj.model,'Flow')
        n = obj.mesh.nNodes;
      elseif isFVTPFABased(obj.model,'Flow')
        n = obj.mesh.nCells;
      end
      obj.state.data.pressure = zeros(n,1);
    end

    function updateState(obj,dSol)
      if nargin > 1
        ents = obj.dofm.getActiveEnts(obj.getField());
        obj.state.data.pressure(ents) = obj.state.data.pressure(ents) + dSol(obj.dofm.getDoF(obj.getField()));
      end
    end

    function var = getState(obj,varargin)
      % input: state structure
      % output: current primary variable
      if isempty(varargin)
        var = obj.state.data.pressure;
      else
        stateIn = varargin{1};
        var = stateIn.data.pressure;
      end
    end

    function setState(obj,id,vals)
      % set values of the primary variable
      obj.state.data.pressure(id) = vals;
    end

    function states = finalizeState(obj,states,t)
      % Compute the posprocessing variables for the module.
      states.potential = computePotential(obj,states.pressure);
      states.head = computePiezHead(obj,states.pressure);
      mob = (1/obj.material.getFluid().getDynViscosity());
      states.flux = computeFlux(obj,states.potential,mob,t);
      states.perm = printPermeab(obj);
      % states.mass = checkMassCons(obj,mob,potential);
    end

    function [cellData,pointData] = printState(obj,sOld,sNew,t)
      % append state variable to output structure
      outPrint = [];
      switch nargin
        case 2
          outPrint.pressure = sOld.data.pressure;
          t = sOld.t;
        case 4
          % linearly interpolate state variables containing print time
          fac = (t - sOld.t)/(sNew.t - sOld.t);
          outPrint.pressure = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);
        otherwise
          error('Wrong number of input arguments');
      end
      outPrint = finalizeState(obj,outPrint,t);
      [cellData,pointData] = SinglePhaseFlow.buildPrintStruct(obj.model,outPrint);
    end

    function computeMat(obj,~,dt)
      % recompute elementary matrices only if the model is linear
      if ~isLinear(obj) || isempty(obj.J)
        if obj.model.isFEMBased('Flow')
          computeMatFEM(obj);
        elseif obj.model.isFVTPFABased('Flow')
          % lw = computeMobility(obj);
          % computeStiffMatFV(obj,lw);
          mu = obj.material.getFluid().getDynViscosity();
          computeStiffMatFV(obj,1/mu);
          computeCapMatFV(obj);
        end
      end
      if obj.simParams.isTimeDependent
        obj.J = obj.simParams.theta*obj.H + obj.P/dt;
      else
        obj.J = obj.H;
      end
    end

    function computeMatFEM(obj)
      % dealing with input params
      subCells = obj.dofm.getFieldCells(obj.getField());
      nEntries = sum(obj.mesh.cellNumVerts(subCells).^2);

      [iiVec,jjVec,HVec,PVec] = deal(zeros(nEntries,1));

      % Get the fluid compressibility
      beta = obj.material.getFluid().getFluidCompressibility();

      % Get the fluid dynamic viscosity
      mu = obj.material.getFluid().getDynViscosity();

      l1 = 0;
      for el = subCells'
        permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
        poro = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPorosity();
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
        dof = obj.dofm.getLocalDoF(nodes,obj.fldId);
        [jjLoc,iiLoc] = meshgrid(dof,dof);
        iiVec(l1+1:l1+s1) = iiLoc(:);
        jjVec(l1+1:l1+s1) = jjLoc(:);
        HVec(l1+1:l1+s1) = HLoc(:);
        PVec(l1+1:l1+s1) = PLoc(:);
        l1 = l1 + s1;
      end
      % renumber indices according to active nodes
      nDoF = obj.dofm.getNumDoF(obj.getField());
      % Assemble H and P matrices defined as new fields of
      obj.H = sparse(iiVec, jjVec, HVec, nDoF, nDoF);
      obj.P = sparse(iiVec, jjVec, PVec, nDoF, nDoF);
    end

    function computeStiffMatFV(obj,lw)
      % Inspired by MRST
      % subCells =
      subCells = obj.dofm.getFieldCells(obj.getField());
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
      nDoF = obj.dofm.getNumDoF(obj.getField());
      obj.H = sparse([neigh1; neigh2; (1:nSubCells)'],...
        [neigh2; neigh1; (1:nSubCells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], nDoF, nDoF);
    end

    function computeCapMatFV(obj,varargin)
      subCells = obj.dofm.getFieldCells(obj.getField());
      nSubCells = length(subCells);
      poroMat = zeros(nSubCells,1);
      alphaMat = zeros(nSubCells,1);
      beta = obj.material.getFluid().getFluidCompressibility();
      for m = 1:obj.mesh.nCellTag
        if ~ismember(m,obj.dofm.getFieldCellTags({obj.getField(),'Poromechanics'}))
          % compute alpha only if there's no coupling in the
          % subdomain
          alphaMat(m) = obj.material.getMaterial(m).ConstLaw.getRockCompressibility();
        end
        poroMat(m) = obj.material.getMaterial(m).PorousRock.getPorosity();
      end
      % (alpha+poro*beta)
      PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));
      if ~isempty(varargin)
        % variably saturated flow model
        PVal = PVal.*varargin{1} + poroMat(obj.mesh.cellTag(subCells)).*varargin{2};
      end
      PVal = PVal.*obj.mesh.cellVolume(subCells);
      nDoF = obj.dofm.getNumDoF(obj.getField());
      [~,~,dof] = unique(subCells);
      obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
    end

    function computeRhs(obj,stateOld,dt)
      % Compute the residual of the flow problem
      lw = 1/obj.material.getFluid().getDynViscosity();
      ents = obj.dofm.getActiveEnts(obj.getField());
      if ~obj.simParams.isTimeDependent
        obj.rhs = obj.H*obj.state.data.pressure(ents);
      else
        theta = obj.simParams.theta;
        rhsStiff = theta*obj.H*obj.state.data.pressure(ents) + (1-theta)*obj.H*stateOld.data.pressure(ents);
        rhsCap = (obj.P/dt)*(obj.state.data.pressure(ents) - stateOld.data.pressure(ents));
        obj.rhs = rhsStiff + rhsCap;
      end
      gamma = obj.material.getFluid().getFluidSpecWeight();
      %adding gravity rhs contribute
      if gamma > 0
        if isFEMBased(obj.model,'Flow')
          obj.rhs = obj.rhs + obj.rhsGrav;
        elseif isFVTPFABased(obj.model,'Flow')
          obj.rhs = obj.rhs + finalizeRHSGravTerm(obj,lw);
        end
      end
    end

    function computeRHSGravTerm(obj)
      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity'
      rhsTmp = zeros(obj.dofm.getNumDoF(obj.getField()),1);
      % rhsTmp = zeros(obj.dofm.getNumDoF(obj.getField()),1);
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0
        subCells = obj.dofm.getFieldCells(obj.getField());
        if isFEMBased(obj.model,'Flow')
          for el = subCells'
            % Get the material permeability
            permMat = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
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
            %
            entsId = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
            rhsTmp(entsId) = rhsTmp(entsId) + rhsLoc;
          end
          obj.rhsGrav = rhsTmp(obj.dofm.getActiveEnts(obj.getField()));
        elseif isFVTPFABased(obj.model,'Flow')
          neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
          zVec = obj.mesh.cellCentroid(:,3);
          zNeigh = zVec(neigh);
          obj.rhsGrav = gamma*obj.trans(obj.isIntFaces).*(zNeigh(:,1) - zNeigh(:,2));
        end
      end
      % remove inactive components of rhs vector
    end

    function gTerm = finalizeRHSGravTerm(obj,lw)
      nCells = obj.dofm.getNumDoF(obj.getField());
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
        -lw.*obj.rhsGrav],[nCells,1]);
      gTerm = gTerm(obj.dofm.getActiveEnts(obj.getField()));
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
          if isFVTPFABased(obj.model,'Flow')
            faceID = obj.bcs.getEntities(id);
            ents = sum(obj.faces.faceNeighbors(faceID,:),2);

            % [ents,~,ind] = unique(ents);
            % % % [faceID, faceOrder] = sort(obj.bcs.getEntities(id));
            % % % ents = sum(obj.faces.faceNeighbors(faceID,:),2);
            % % % v(faceOrder,1) = obj.bcs.getVals(id,t);
            switch obj.bcs.getType(id)
              case 'Neu'
                vals = vecnorm(obj.faces.faceNormal(faceID,:),2,2).*v;
              case 'Dir'
                gamma = obj.material.getFluid().getFluidSpecWeight();
                mu = obj.material.getFluid().getDynViscosity();
                tr = obj.getFaceTransmissibilities(faceID);

                % q = 1/mu*tr.*((obj.state.data.pressure(ents) - v)...
                %    + gamma*(obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3)));
                % vals = [1/mu*tr,accumarray(ind,q)]; % {JacobianVal,rhsVal]

                dirJ = 1/mu*tr;
                % % press = obj.state.data.pressure(ents) - v;
                % % gravT =  gamma*(obj.mesh.cellCentroid(ents,3) ...
                % %   - obj.faces.faceCentroid(faceID,3));
                potential = (obj.state.data.pressure(ents) - v) ...
                  + gamma*(obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3));
                q = dirJ.*potential;
                vals = [dirJ,q];

              case 'Spg'
                gamma = obj.material.getFluid().getFluidSpecWeight();
                assert(gamma>0.,'To impose Seepage boundary condition is necessary the fluid specify weight be bigger than zero!');

                % Datum = max(obj.mesh.coordinates);
                % zbc = Datum(3)-obj.faces.faceCentroid(faceID,3);
                % href = Datum(3)-bc.getVals(id,t);
                % href = Datum(3)-href(1);
                % v = gamma*(zbc-href);
                % zbc = obj.faces.faceCentroid(faceID,3);
                % href = bc.getVals(id,t);
                % v = gamma*(href(1)-zbc);
                zbc = obj.faces.faceCentroid(faceID,3);
                href = v(1);
                v = gamma*(href-zbc);

                v(v<=0)=0.;
                mu = obj.material.getFluid().getDynViscosity();
                tr = obj.getFaceTransmissibilities(faceID);
                q = 1/mu*tr.*((obj.state.data.pressure(ents) - v)...
                  + gamma*(obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3)));
                vals = [1/mu*tr,q];
            end
          elseif isFEMBased(obj.model,'Flow')
            ents = obj.bcs.getLoadedEntities(id);
            entitiesInfl = obj.bcs.getEntitiesInfluence(id);
            vals = entitiesInfl*v;
          end
        case 'VolumeForce'
          v = obj.bcs.getVals(id,t);
          if isFVTPFABased(obj.model,'Flow')
            ents = obj.bcs.getEntities(id);
            vals = v.*obj.mesh.cellVolume(ents);
          elseif isFEMBased(obj.model,'Flow')
            ents = obj.bcs.getLoadedEntities(id);
            entitiesInfl = obj.bcs.getEntitiesInfluence(id);
            vals = entitiesInfl*v;
          end
      end
      % get local dof numbering
      dof = obj.dofm.getLocalDoF(ents,obj.fldId);
    end

    function applyDirVal(obj,dof,vals)
      if isFVTPFABased(obj.model,'Flow')
        % Dirichlet BCs cannot be directly applied to the solution
        % vector
        return
      end
      obj.state.data.pressure(dof) = vals;
    end

    function applyDirBC(obj,~,ents,vals)
      % apply Dirichlet BCs
      % ents: id of constrained faces without any dof mapping applied
      % vals(:,1): Jacobian BC contrib vals(:,2): rhs BC contrib
      if isFVTPFABased(obj.model,'Flow')
        % BCs imposition for finite volumes
        assert(size(vals,2)==2,'Invalid matrix size for BC values');
        nDoF = obj.dofm.getNumDoF(obj.getField());
        obj.J(nDoF*(ents-1) + ents) = obj.J(nDoF*(ents-1) + ents) + vals(:,1);
        obj.rhs(ents) = obj.rhs(ents) + vals(:,2);
      else
        % strong nodal BCs imposition
        applyDirBC@SinglePhysics(obj,obj.getField(),ents)
      end
    end

    function out = isLinear(obj)
      out = true;
    end

    function trans = getFaceTransmissibilities(obj,faceID)
      trans = obj.trans(faceID);
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
        KMat(i,:) = obj.material.getMaterial(i).PorousRock.getPermVector();
      end
      hT = zeros(length(hf2Cell),1);
      for k=1:length(r)
        hT = hT + L(:,r(k)) .* KMat(obj.mesh.cellTag(hf2Cell),k) .* N(:,c(k));
      end
      hT = hT./sum(L.*L,2);
      %       mu = obj.material.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
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
        KMat(i,:) = obj.material.getMaterial(i).PorousRock.getPermVector();
      end
      hT = zeros(length(hf2Cell),1);
      for k=1:length(r)
        hT = hT + L(:,r(k)) .* KMat(obj.mesh.cellTag(hf2Cell),k) .* N(:,c(k));
      end
      trans = hT./sum(L.*L,2);
    end

    function alpha = getRockCompressibility(obj,el)
      if ismember(obj.mesh.cellTag(el),getFieldCellTags(obj.dofm,{obj.getField(),'Poromechanics'}))
        alpha = 0; %this term is not needed in coupled formulation
      else
        if isfield(obj.material.getMaterial(obj.mesh.cellTag(el)),"ConstLaw")
          %solid skeleton contribution to storage term as oedometric compressibility .
          alpha = obj.material.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
        else
          alpha = 0;
        end
      end
    end

    function pzhead = computePiezHead(obj,pressure)
      % COMPUTEPIEZHEAD - compute the piezometric head for the cell or element.
      if isFEMBased(obj.model,'Flow')
        zbc=obj.mesh.coordinates(:,3);
      elseif isFVTPFABased(obj.model,'Flow')
        zbc=obj.mesh.cellCentroid(:,3);
      end
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0.
        pzhead = zbc+pressure/gamma;
      else
        pzhead = zeros(length(zbc),1);
      end
    end

    function potential = computePotential(obj,pressure)
      % COMPUTEFLUX - compute the potential for the cell or element.
      potential = pressure;
      gamma = obj.material.getFluid().getFluidSpecWeight();
      if gamma > 0
        if isFEMBased(obj.model,'Flow')
          potential = potential + gamma*obj.mesh.coordinates(:,3);
        elseif isFVTPFABased(obj.model,'Flow')
          potential = potential + gamma*obj.mesh.cellCentroid(:,3);
        end
      end
    end

    function perm = printPermeab(obj)
      %printPropState - print the potential for the cell or element.
      perm = zeros(obj.mesh.nCells,6);
      for el=1:obj.mesh.nCells
        ktmp = obj.material.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
        perm(el,1)=ktmp(1,1);
        perm(el,2)=ktmp(2,2);
        perm(el,3)=ktmp(3,3);
        perm(el,4)=ktmp(1,2);
        perm(el,5)=ktmp(2,3);
        perm(el,6)=ktmp(1,3);
      end
    end

    function flux = computeFlux(obj,pot,mob,t)
      %COMPUTEFLUX - compute the flux at the faces, than accumulate
      %the value at the nodes (The contribution of the boundary is done
      % in another function).
      flux = zeros(obj.mesh.nNodes,3);
      if isFEMBased(obj.model,'Flow') & false
        % TODO - This part is still need some work.
      elseif isFVTPFABased(obj.model,'Flow')
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
          if strcmp(obj.model.findPhysicsFromID(obj.model.findIDPhysics(obj.bcs.getPhysics(bc))),'Flow')
            [~,vals] = getBC(obj,bc,t);
            [ents, ~] = sort(obj.bcs.getEntities(bc));
            switch obj.bcs.getCond(bc)
              case {'NodeBC','ElementBC'}
              case 'SurfBC'
                if strcmp(obj.bcs.getType(bc),'Dir') || strcmp(obj.bcs.getType(bc),'Spg')
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
    end

      % function mass = checkMassCons(obj,mob,pot)
      %   %CHECKMASSCONS - check the mass conservation in all elements.
      %   mass = zeros(obj.mesh.nCells,1);
      %   if isFVTPFABased(obj.model,'Flow')
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
        mu = obj.material.getFluid().getDynViscosity();
        if mu==0
          lwkpt = 1;
        else
          lwkpt = 1/mu;
        end
        dlwkpt = 0;
      end
    end


    methods (Static)
      function [cellStr,pointStr] = buildPrintStruct(mod,state)
        if isFEMBased(mod,'Flow')
          cellStr = repmat(struct('name', 1, 'data', 1), 1, 1);
          cellStr(1).name = 'permeability';
          cellStr(1).data = state.perm;

          pointStr = repmat(struct('name', 1, 'data', 1), 2, 1);
          pointStr(1).name = 'pressure';
          pointStr(1).data = state.pressure;
          pointStr(2).name = 'potential';
          pointStr(2).data = state.potential;
          pointStr(3).name = 'piezometric head';
          pointStr(3).data = state.head;
        elseif isFVTPFABased(mod,'Flow')
          pointStr = repmat(struct('name', 1, 'data', 1), 1, 1);
          pointStr(1).name = 'flux';
          pointStr(1).data = state.flux;
          % pointStr = [];

          cellStr = repmat(struct('name', 1, 'data', 1), 3, 1);
          cellStr(1).name = 'permeability';
          cellStr(1).data = state.perm;
          cellStr(2).name = 'pressure';
          cellStr(2).data = state.pressure;
          cellStr(3).name = 'potential';
          cellStr(3).data = state.potential;
          cellStr(4).name = 'piezometric head';
          cellStr(4).data = state.head;
        end
      end

      function out = getField()
        out = SinglePhaseFlow.field;
      end
    end
  end

