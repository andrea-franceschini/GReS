classdef SinglePhaseFlow < PhysicsSolver
  %SINGLEPHASEFLOW


  properties
    trans
    isIntFaces
    rhsGrav         % gravity contribution to rhs
    H
    P
    scheme = "FiniteVolumesTPFA"        % or FEM
  end

  properties (Access = private)
    fieldId
  end

  methods (Access = public)
    function obj = SinglePhaseFlow(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,solverInput)

      nTags = obj.mesh.nCellTag;

      if ~isempty(solverInput)
        targetRegions = getXMLData(solverInput,1:nTags,"targetRegions");
        obj.scheme = getXMLData(solverInput,obj.scheme,"scheme");
      else
        targetRegions = 1:nTags;
      end


      switch obj.scheme
        case "FiniteVolumesTPFA"
          obj.dofm.registerVariable(obj.getField(),entityField.cell,1,targetRegions);
        case "FEM"
          obj.dofm.registerVariable(obj.getField(),entityField.node,1,targetRegions);
        otherwise
          error("Scheme %s for class %s is not a valid GReS scheme",...
            obj.scheme,class(obj));
      end

      obj.fieldId = obj.dofm.getVariableId(obj.getField());

      n = getNumbDoF(obj.dofm,obj.fieldId);

      % initialize the state object with a pressure field
      obj.getState().data.(obj.getField()) = zeros(n,1);

      if strcmp(obj.scheme,'FiniteVolumesTPFA')

        linkBoundSurf2TPFAFace(obj);

        obj.computeTrans();
        %get cells with active flow model
        flowCells = obj.dofm.getActiveEntities(obj.fieldId);
        % Find internal faces (i.e. shared by two active flow cells)
        obj.isIntFaces = all(ismember(obj.faces.faceNeighbors, flowCells), 2);
      end

      computeRHSGravTerm(obj);

    end


    function assembleSystem(obj,dt)

      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);

      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);

    end
    %


    function updateState(obj,dSol)
      if nargin > 1
        ents = obj.dofm.getActiveEntities(obj.fieldId);
        state = getState(obj);
        state.data.pressure(ents) = state.data.pressure(ents) + dSol(obj.dofm.getDoF(obj.fieldId));
      end
    end

    function advanceState(obj)
      % does nothing for now, but needed to override the abstract
      % physicsSolver method
    end


    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      states.potential = computePotential(obj,p);
      states.head = computePiezHead(obj,p);
      mob = (1/obj.materials.getFluid().getDynViscosity());
      states.flux = computeFlux(obj,p,mob,t);
      states.perm = printPermeab(obj);
      states.pressure = p;
      % states.mass = checkMassCons(obj,mob,potential);
    end

    function [cellData,pointData] = printVTK(obj,t)
      % append state variable to output structure

      sOld = getStateOld(obj);
      sNew = getState(obj);

      fac = (t - sOld.t)/(sNew.t - sOld.t);
      p = sNew.data.pressure*fac + sOld.data.pressure*(1-fac);
      outPrint = finalizeState(obj,p,t);
      [cellData,pointData] = SinglePhaseFlow.buildPrintStruct(outPrint);

    end


    function J = computeMat(obj,dt)

      % recompute elementary matrices only if the model is linear
      id = getVariableId(obj.dofm,obj.fieldId);
      if ~isLinear(obj) || isempty(obj.J{id,id})
        if isFEM(obj)
          computeMatFEM(obj);
        elseif isTFPA(obj)
          mu = obj.materials.getFluid().getDynViscosity();
          computeStiffMatFV(obj,1/mu);
          computeCapMatFV(obj);
        end
      end

      if obj.simparams.isTimeDependent
        J = obj.simparams.theta*obj.H + obj.P/dt;
      else
        J = obj.H;
      end
    end

    function computeMatFEM(obj)
      % dealing with input params
      subCells = obj.dofm.getTargetRegions(obj.fieldId);
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
        dof = obj.dofm.getLocalDoF(nodes,obj.fieldId);
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

    function computeStiffMatFV(obj,lw)
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
      nDoF = obj.dofm.getNumDoF(obj.fieldId);
      obj.H = sparse([neigh1; neigh2; (1:nSubCells)'],...
        [neigh2; neigh1; (1:nSubCells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], nDoF, nDoF);
    end

    function computeCapMatFV(obj,varargin)
      subCells = obj.dofm.getFieldCells(obj.fieldId);
      nSubCells = length(subCells);
      poroMat = zeros(nSubCells,1);
      alphaMat = zeros(nSubCells,1);
      beta = obj.materials.getFluid().getFluidCompressibility();
      for m = 1:obj.mesh.nCellTag
        if ~ismember(m,obj.dofm.getTargetRegions([obj.fieldId,'displacements']))
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
      p = getState(obj,obj.fieldId);
      pOld = getStateOld(obj,obj.fieldId);

      lw = 1/obj.materials.getFluid().getDynViscosity();
      ents = obj.dofm.getActiveEnts(obj.fieldId);

      if ~obj.simparams.isTimeDependent
        rhs = obj.H*p(ents);
      else
        theta = obj.simparams.theta;
        rhsStiff = theta*obj.H*p(ents) + (1-theta)*obj.H*pOld(ents);
        rhsCap = (obj.P/dt)*(p(ents) - pOld(ents));
        rhs = rhsStiff + rhsCap;
      end

      %adding gravity rhs contribute
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma > 0
        if isFEM(obj)
          rhs = rhs + obj.rhsGrav;
        elseif isTPFA(obj)
          rhs = rhs + finalizeRHSGravTerm(obj,lw);
        end
      end

    end

    % TO DO: update to new FEM logic
    function computeRHSGravTerm(obj)
      % Compute the gravity contribution
      % Get the fluid specific weight and viscosity'
      rhsTmp = zeros(obj.dofm.getNumbDoF(obj.fieldId),1);
      % rhsTmp = zeros(obj.dofm.getNumDoF(obj.fieldId),1);
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma > 0
        subCells = obj.dofm.getFieldCells(obj.fieldId);
        if isFEM(obj)
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
            %
            entsId = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
            rhsTmp(entsId) = rhsTmp(entsId) + rhsLoc;
          end
          obj.rhsGrav = rhsTmp(obj.dofm.getActiveEnts(obj.fieldId));
        elseif isTPFA(obj)
          neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
          zVec = obj.mesh.cellCentroid(:,3);
          zNeigh = zVec(neigh);
          obj.rhsGrav = gamma*obj.trans(obj.isIntFaces).*(zNeigh(:,1) - zNeigh(:,2));
        end
      end
      % remove inactive components of rhs vector
    end

    function gTerm = finalizeRHSGravTerm(obj,lw)
      nCells = obj.dofm.getNumDoF(obj.fieldId);
      neigh = obj.faces.faceNeighbors(obj.isIntFaces,:);
      gTerm = accumarray(neigh(:),[lw.*obj.rhsGrav; ...
        -lw.*obj.rhsGrav],[nCells,1]);
      gTerm = gTerm(obj.dofm.getActiveEnts(obj.fieldId));
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
          if isTPFA(obj)
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
                tr = obj.getFaceTransmissibilities(faceID);

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
                tr = obj.getFaceTransmissibilities(faceID);
                dz = obj.mesh.cellCentroid(ents,3) - obj.faces.faceCentroid(faceID,3);
                q = 1/mu*tr.*(p(ents) - v) + gamma*dz;
                vals = [1/mu*tr,q];
            end

          elseif isFEM(obj)
            ents = obj.bcs.getLoadedEntities(id);
            entitiesInfl = obj.bcs.getEntitiesInfluence(id);
            vals = entitiesInfl*v;
          end

        case 'VolumeForce'
          v = obj.bcs.getVals(id,t);
          if isTPFA(obj)
            ents = obj.bcs.getEntities(id);
            vals = v.*obj.mesh.cellVolume(ents);
          elseif isFEM(obj)
            ents = obj.bcs.getLoadedEntities(id);
            entitiesInfl = obj.bcs.getEntitiesInfluence(id);
            vals = entitiesInfl*v;
          end
      end

      dof = obj.dofm.getLocalDoF(obj.fieldId,ents);

    end

    function applyDirVal(obj,bcId,t)
      if isTPFA(obj)
        % Dirichlet BCs cannot be directly applied to the solution
        % vector
        return
      end

      [bcDofs,bcVals] = getBC(obj,bcId,t);

      state = getState(obj);
      state.data.pressure(bcDofs) = bcVals;
    end

    function applyBC(obj,t,bcId,bcVar)

      [bcDofs,bcVals] = getBC(obj,bcId,t);

      % Base application of a Boundary condition
      bcType = obj.bcs.getType(bcId);

      switch bcType
        case {'Dirichlet','Seepage'}
          applyDirBC(obj,bcDofs,bcVals,bcVar);
        case 'Neumann'
          applyNeuBC(obj,bcDofs,bcVar);
        otherwise
          error("Error in %s: Boundary condition type '%s' is not " + ...
            "available in applyBC()",class(obj),bcType)
      end
    end


    function [cellData,pointData] = writeVTK(obj,t)

      % append state variable to output structure
      sOld = getStateOld(obj);
      sNew = getState(obj);

      if isempty(sOld)
        p = sNew.data.pressure;
      else
        fac = (t - sOld.t)/(sNew.t - sOld.t);
        p = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);
      end

      outPrint = finalizeState(obj,p,t);
      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function writeMatFile(obj,t)

        pOld = getStateOld(obj,obj.getField());
        pCurr = getState(obj,obj.getField());
        fac = (t - getState(obj).t)/(getState(obj).t - getStateOld(obj).t);
        obj.outstate.results(tID+1).pressure = pCurr*fac+pOld*(1-fac);

    end




    function applyDirBC(obj,bcDofs,bcVals,bcVar)
      % apply Dirichlet BCs
      % ents: id of constrained faces without any dof mapping applied
      % vals(:,1): Jacobian BC contrib vals(:,2): rhs BC contrib
      if isTPFA(obj)
        % BCs imposition for finite volumes - boundary flux
        assert(size(bcVals,2)==2,'Invalid matrix size for BC values');
        nDoF = obj.dofm.getNumDoF(obj.fieldId);
        obj.J(nDoF*(bcDofs-1) + bcDofs) = obj.J(nDoF*(bcDofs-1) + bcDofs) + bcVals(:,1);
        obj.rhs(bcDofs) = obj.rhs(bcDofs) + bcVals(:,2);
      else
        % FEM - strong BCs imposition
        applyDirBC@PhysicsSolver(obj,bcDofs,bcVar)
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

    function alpha = getRockCompressibility(obj,el)
      if ismember(obj.mesh.cellTag(el),getTargetRegions(obj.dofm,["pressure","displacements"]));
        alpha = 0; %this term is not needed in a coupled formulation
      else
        if isfield(obj.materials.getMaterial(obj.mesh.cellTag(el)),"ConstLaw")
          %solid skeleton contribution to storage term as oedometric compressibility .
          alpha = obj.materials.getMaterial(obj.mesh.cellTag(el)).ConstLaw.getRockCompressibility();
        else
          alpha = 0;
        end
      end
    end

    function pzhead = computePiezHead(obj,pressure)
      % COMPUTEPIEZHEAD - compute the piezometric head for the cell or element.
      if isFEM(obj)
        zbc = obj.mesh.coordinates(:,3);
      elseif isTPFA(obj)
        zbc = obj.mesh.cellCentroid(:,3);
      end
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma > 0.
        pzhead = zbc+pressure/gamma;
      else
        pzhead = zeros(length(zbc),1);
      end
    end

    function potential = computePotential(obj,pressure)
      % COMPUTEFLUX - compute the potential for the cell or element.
      potential = pressure;
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma > 0
        if isFEM(obj)
          potential = potential + gamma*obj.mesh.coordinates(:,3);
        elseif isTPFA(obj)
          potential = potential + gamma*obj.mesh.cellCentroid(:,3);
        end
      end
    end

    function perm = printPermeab(obj)
      %printPropState - print the potential for the cell or element.
      perm = zeros(obj.mesh.nCells,6);
      for el=1:obj.mesh.nCells
        ktmp = obj.materials.getMaterial(obj.mesh.cellTag(el)).PorousRock.getPermMatrix();
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
      if isFEM(obj) & false
        % TODO - This part still need some work.
      elseif isTPFA(obj)
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
    end

    function out = isFEM(obj)
      out = false;
      if strcmp(obj.scheme,"FEM")
        out = true;
      end
    end

    function out = isTPFA(obj)
      out = false;
      if strcmp(obj.scheme,"FiniteVolumesTPFA")
        out = true;
      end
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

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      if isFEM(obj)
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
      elseif isTPFA(obj)
        pointStr = repmat(struct('name', 1, 'data', 1), 1, 1);
        pointStr(1).name = 'flux';
        pointStr(1).data = state.flux;

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

  end



  methods (Static)

    function out = getField()
      out = "pressure";

    end

  end

end
