classdef SolidMechanicsContact < MeshTying

  % Mortar contact with piecewise constant traction multipliers.
  % isSmooth = true: fixed active set in Newton, updated in an outer loop.
  % isSmooth = false: contact branches selected inside each Newton assembly.

  properties
    phi               % friction angle in deg
    cohesion          % cohesion
    contactHelper
    activeSet
    NLIter = 0
    stickNodes     % boundary nodes where contact state should stay stick
    forceStick     % flag to enforce interface to stay stick
    contactAugmentation
    isSmooth = true % outer active-set loop; false selects semi-smooth Newton
  end

  methods

    function obj = SolidMechanicsContact(id,domains,inputStruct)

      obj@MeshTying(id,domains,inputStruct);

      if obj.multiplierLocation ~= entityField.surface
        error("Interface Solver %s is not implemented for multipliers" + ...
          "located at %s. The only available entityField is %s",...
          class(obj),obj.multiplierLocation,entityField.surface)
      end

    end

    function registerInterface(obj,varargin)

      input = varargin{1};

      input = readInput(struct('Coulomb',[],'ActiveSet',missing,'forceStick',0,...
        'stabilizationScale',1.0,'augmentationParameter',1.0,...
        'augmentationNormal',1.0,'augmentationTangential',1.0,'isSmooth',1),input);

      params = readInput(struct('cohesion',[],'frictionAngle',[]),input.Coulomb);

      obj.stabilizationScale = input.stabilizationScale;

      obj.contactAugmentation(1) = input.augmentationNormal;
      obj.contactAugmentation(2) = input.augmentationTangential;

      obj.isSmooth = logical(input.isSmooth);

      obj.forceStick = logical(input.forceStick);

      obj.cohesion = params.cohesion;
      obj.phi = params.frictionAngle;

      nDofsInterface = getNumbDoF(obj);

      s = getState(obj);

      s.traction = zeros(nDofsInterface,1);
      s.deltaTraction = zeros(nDofsInterface,1);

      % Raw mortar gap in the local contact frame
      s.gap = zeros(nDofsInterface,1);
      s.normalGap = zeros(round(1/3*nDofsInterface),1);

      % total tangential gap
      s.tangentialGap = zeros(round(2/3*nDofsInterface),1);

      % gap variation at the current time step
      s.tangentialSlip = zeros(round(2/3*nDofsInterface),1);
      setState(obj,s);

      N = obj.grids(MortarSide.slave).surfaces.num;
      initializeActiveSet(obj,N,input.ActiveSet);
      
    end

    function updateState(obj,du)

      % traction update
      actMult = getMultiplierDoF(obj);

      state = getState(obj);
      stateOld = getStateOld(obj);
      state.traction(actMult) = state.traction(actMult) + du(1:obj.nMult);
      state.deltaTraction = state.traction - stateOld.traction;
      obj.NLIter = obj.NLIter + 1;

      setState(obj,state);

      % The smooth stick trial must retain inadmissible reaction tractions
      % so that the outer loop can detect opening/sliding.
      % if ~obj.isSmooth
      %   applyContactReturnMap(obj);
      % end

      % update gap
      computeGap(obj);

      if gresLog().getVerbosity >= 2
        nStick = sum(obj.activeSet.curr == ContactMode.stick);
        nSlip = sum(obj.activeSet.curr == ContactMode.slip | ...
          obj.activeSet.curr == ContactMode.newSlip);
        nOpen = sum(obj.activeSet.curr == ContactMode.open);

        fprintf('%s: active set ',class(obj));
        fprintf('(NLIter %i): stick = %i, slip = %i, open = %i\n', ...
          obj.NLIter,nStick,nSlip,nOpen);
      end

    end

    function assembleConstraint(obj)

      % reset the jacobian blocks
      obj.setJmu(MortarSide.slave, []);
      obj.setJmu(MortarSide.master, []);
      obj.setJum(MortarSide.slave, []);
      obj.setJum(MortarSide.master, []);

      if isempty(obj.D)
        computeConstraintMatrices(obj);
      end

      computeContactMatricesAndRhs(obj);

      % Semi-smooth stabilization is already differentiated inside the law.
      if ~obj.isSmooth
        return
      end

      % get stabilization matrix depending on the current active set
      [H,rhsStab] = getStabilizationMatrixAndRhs(obj);

      obj.Jconstraint = obj.Jconstraint - H;
      obj.rhsConstraint = obj.rhsConstraint + rhsStab;

      if gresLog().getVerbosity > 3
        % print rhs terms for each fracture state for debug purposes
        dof_stick = DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.stick),3);
        dof_slip = [DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.slip),3); ...
                    DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.newSlip),3)];
        dof_open = DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.open),3);
        fprintf('Rhs norm for stabilization: %4.3e \n', norm(rhsStab));
        fprintf('Rhs norm for stick dofs: %4.3e \n', norm(obj.rhsConstraint(dof_stick)))
        fprintf('Rhs norm for slip dofs: %4.3e \n', norm(obj.rhsConstraint(dof_slip)))
        fprintf('Rhs norm for open dofs: %4.3e \n', norm(obj.rhsConstraint(dof_open)))
      end

    end

    function applyContactReturnMap(obj)

      % Semi-smooth traction postprocessing, as in the Augmented class.

      state = getState(obj);
      stateOld = getStateOld(obj);
      tanPhi = tan(deg2rad(obj.phi));

      for is = 1:numel(obj.activeSet.curr)

        % Forced-stick constraints may support reactions outside Coulomb.
        if isForceStickElement(obj,is)
          continue
        end
        id = DoFManager.dofExpand(is,3);
        tTrial = state.traction(id);

        % normal projection
        tN = min(tTrial(1),0.0);

        % Coulomb bound
        tauLim = max(obj.cohesion - tanPhi*tN,0.0);

        % tangential projection
        tTTrial = tTrial(2:3);
        tTNorm = norm(tTTrial);

        if tTNorm > tauLim && tTNorm > 0
          tT = tauLim * tTTrial / tTNorm;
        else
          tT = tTTrial;
        end

        state.traction(id) = [tN; tT];

      end

      state.deltaTraction = state.traction - stateOld.traction;
      setState(obj,state);

    end

    function hasConfigurationChanged = updateConfiguration(obj)
      % Only the smooth strategy asks the driver for another outer solve.
      if ~obj.isSmooth
        hasConfigurationChanged = false;
        obj.NLIter = 0;
        return
      end
      hasConfigurationChanged = updateOuterActiveSet(obj);
    end

    function initialize(obj)

      initialize@InterfaceSolver(obj);

      % initial traction from cell stress
      tIni = computeInitialTraction(obj);

      addInitialTraction(obj,tIni);

      setStateOld(obj,getState(obj,"traction"),"traction");

      setStickNodes(obj);
      enforceForcedStick(obj);

    end

    function timeStepSetup(obj)
      obj.NLIter = 0;
      if obj.isSmooth
        % Start the outer active-set iteration with a stick trial.
        isActive = obj.activeSet.curr ~= ContactMode.open;
        obj.activeSet.curr(isActive) = ContactMode.stick;
      end
      enforceForcedStick(obj);
    end

    function out = isLinear(obj)
      % The fixed all-stick contact equations are linear. Semi-smooth
      % classification can change during Newton even if currently all stick.
      out = obj.isSmooth && all(obj.activeSet.curr == ContactMode.stick);
    end

    function addInitialTraction(obj,tIni)
      % add a traction on the fault
      t = getState(obj,"traction");
      t = t + tIni;
      setState(obj,t,"traction");
      setStateInit(obj,t,"traction");

    end

    function trac = computeInitialTraction(obj)
      % initialize traction for cell stress (average)
      sl = MortarSide.slave;
      avgStress = obj.domains(MortarSide.slave).getState("avgStress");
      surf = obj.grids(sl).surfaces;
      faces = obj.domains(sl).grid.faces;
      normals = surf.normal;
      faceIds = surf.faceId;
      cellIds = faces.neighbors(faceIds,1);
      sigma = zeros(3);
      idx = [1;6;5;6;2;4;5;4;3];
      trac = zeros(getNumbDoF(obj),1);
      for i = 1:numel(cellIds)
        cellId = cellIds(i);
        sigma(:) = avgStress(cellId,idx);
        n = normals(i,:);
        tDof = getMultiplierDoF(obj,i);
        t = sigma*n';   % global
        R = getRotationMatrix(obj,sl,i);
        trac(tDof) = R'*t;
      end
    end

    function advanceState(obj)

      advanceState@InterfaceSolver(obj);
      state = getState(obj);
      state.deltaTraction(:) = 0;
      setState(obj,state);
      obj.activeSet.prev = obj.activeSet.curr;
      obj.NLIter = 0;

      % reset the counter for changed states
      obj.activeSet.stateChange(:) = 0;

    end

    function isReset = resetConfiguration(obj)

      toReset = obj.activeSet.curr(:) ~= ContactMode.open;
      obj.activeSet.curr(toReset) = ContactMode.stick;

      enforceForcedStick(obj);
      isReset = true;
    end

    function goBackState(obj,dt)

      % reset state to beginning of time step
      goBackState@InterfaceSolver(obj);
      state = getState(obj);
      state.deltaTraction(:) = 0;
      setState(obj,state);

      obj.activeSet.curr = obj.activeSet.prev;
      obj.NLIter = 0;
      if obj.activeSet.resetActiveSet
        resetConfiguration(obj);
      end
      enforceForcedStick(obj);
    end

    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      outTraction = obj.state.interpolate(fac,"traction");
      dT = obj.state.interpolate(fac,"deltaTraction");
      outNormalGap = obj.state.interpolate(fac,"normalGap");
      outTangentialSlip = obj.state.interpolate(fac,"tangentialSlip");
      outTangentialGap = obj.state.interpolate(fac,"tangentialGap");

      outTangentialSlip = (reshape(outTangentialSlip,2,[]))';
      outTangentialGap = (reshape(outTangentialGap,2,[]))';

      outTangentialGapNorm = sqrt(outTangentialGap(:,1).^2 + ...
        outTangentialGap(:,2).^2);

      tT = [outTraction(2:3:end),outTraction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      deltaTrac = [dT(1:3:end), dT(2:3:end), dT(3:3:end)];

      fractureState = double(obj.activeSet.curr);

      pointStr = [];

      entries = {
        'normal_gap',              outNormalGap
        'normal_stress',           outTraction(1:3:end)
        'tangential_traction_1',   outTraction(2:3:end)
        'tangential_traction_2',   outTraction(3:3:end)
        'tangential_traction_norm',norm_tT
        'tangential_slip',         outTangentialSlip
        'tangential_gap',          outTangentialGap
        'tangential_gap_norm',     outTangentialGapNorm
        'fracture_state',          fractureState
        'rotationMatrix',          obj.grids(1).surfaces.rotationMatrices
        'deltaTraction',           deltaTrac
        };

      surfaceStr = cell2struct(entries, {'name','data'}, 2);
    end

    function writeSolution(obj,fac,tID)

      s = obj.state.interpolate(fac);

      tT = [s.traction(2:3:end),s.traction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      obj.outstate.results(tID).tractionVec = s.traction;
      obj.outstate.results(tID).normalGap = s.normalGap;
      obj.outstate.results(tID).slipIncrement = s.tangentialSlip;
      obj.outstate.results(tID).tangentialGap = s.tangentialGap;
      obj.outstate.results(tID).tangentialTractionNorm = norm_tT;

    end

  end

  methods (Access = protected)

    function hasConfigurationChanged = updateOuterActiveSet(obj)
      hasConfigurationChanged = false;
      % Update the fixed active set after the inner Newton solve.

      if obj.forceStick
        return
      end

      obj.NLIter = 0;

      oldActiveSet = obj.activeSet.curr;
      surfSlave = obj.grids(MortarSide.slave).surfaces;

      state = getState(obj);

      for is = 1:numel(obj.activeSet.curr)

        currAS = obj.activeSet.curr(is);

        if isForceStickElement(obj,is)
          obj.activeSet.curr(is) = ContactMode.stick;
          continue
        end

        id = DoFManager.dofExpand(is,3);
        t = state.traction(id);
        limitTraction = abs(obj.cohesion - tan(deg2rad(obj.phi))*t(1));

        % report traction during activeSet update
        gresLog().log(5,['\n Element %i: traction: %1.4e %1.4e %1.4e   ' ...
          'Limit tangential traction: %1.4e \n'],is,t(:), limitTraction)

        obj.activeSet.curr(is) = updateContactState(currAS,t,...
          limitTraction, ...
          state.normalGap(is),...
          obj.activeSet.tol);

      end

      % check if active set changed
      asNew = obj.activeSet.curr;
      asOld = oldActiveSet;

      % Do not update an element that exceeded the maximum number of
      % individual updates
      reset = obj.activeSet.stateChange >= ...
        obj.activeSet.tol.maxStateChange;

      asNew(reset) = asOld(reset);

      obj.activeSet.curr = asNew;
      diffState = asNew - asOld;

      idNewSlipToSlip = all([asOld==2 diffState==1],2);
      diffState(idNewSlipToSlip) = 0;
      hasChangedElem = diffState~=0;

      nomoreStick = diffState > 0;

      obj.activeSet.stateChange(nomoreStick) = ...
        obj.activeSet.stateChange(nomoreStick) + 1;

      hasConfigurationChanged = any(diffState);

      gresLog().log(2,'%s: Active set \n',class(obj));

      if gresLog().getVerbosity > 3
        % report active set changes
        da = asNew - asOld;
        d = da(asOld == 1);
        assert(~any(d==2));       % avoid stick to slip without newSlip
        fprintf('%i elements from stick to new slip \n',sum(d==1));
        fprintf('%i elements from stick to open \n',sum(d==3));
        d = da(asOld==2);
        fprintf('%i elements from new slip to stick \n',sum(d==-1));
        fprintf('%i elements from new slip to slip \n',sum(d==1));
        fprintf('%i elements from new slip to open \n',sum(d==2));
        d = da(asOld==3);
        fprintf('%i elements from slip to stick \n',sum(d==-2));
        fprintf('%i elements from slip to open \n',sum(d==1));
        d = da(asOld==4);
        fprintf('%i elements from open to stick \n',sum(d==-3));
      end

      gresLog().log(2,'Stick dofs: %i    Slip dofs: %i    Open dofs: %i \n',...
        sum(asNew==1), sum(any([asNew==2,asNew==3],2)), sum(asNew==4));

      if hasConfigurationChanged

        % EXCEPTION 1): check if area of fracture changing state is relatively small

        areaChanged = sum(surfSlave.area(hasChangedElem));
        totArea = sum(surfSlave.area);
        if areaChanged/totArea < obj.activeSet.tol.areaChange
          % change the active set, but flag it as nothing changed
          hasConfigurationChanged = false;
          gresLog().log(1,['Active set update suppressed due to small fracture change:' ...
            ' areaChange/areaTot = %3.2e \n'],areaChanged/totArea);
        end

        % EXCEPTION 2): check if changing elements have been looping from
        % stick to slip/open too much times

        if all(obj.activeSet.stateChange(hasChangedElem) > obj.activeSet.tol.maxStateChange)
          hasConfigurationChanged = false;
          gresLog().log(1,['Active set update suppressed due to' ...
            ' unstable behavior detected'])
        end
      end
    end

    function computeGap(obj)
      % compute normal gap and tangential slip (local coordinates)

      state = getState(obj);
      stateOld = getStateOld(obj);

      um = obj.domains(MortarSide.master).getState("displacements");
      us = obj.domains(MortarSide.slave).getState("displacements");

      areaSlave = repelem(obj.getSlaveArea(),3,1);
      areaGap = obj.D*us + obj.M*um;
      state.gap = areaGap./areaSlave;

      if obj.isSmooth
        % Preserve the original outer active-set strategy.
        [~,rhsStab] = getStabilizationMatrixAndRhs(obj);
        stabGap = (areaGap + rhsStab)./areaSlave;
        stabSlip = (state.gap-stateOld.gap) + rhsStab./areaSlave;
        stabSlip(1:3:end) = [];
        state.tangentialSlip = stabSlip;
        state.normalGap = stabGap(1:3:end);
        state.tangentialGap = stateOld.tangentialGap + stabSlip;
      else
        % Physical outputs contain only the geometric jump, never H*t.
        rawSlip = state.gap-stateOld.gap;
        rawSlip(1:3:end) = [];
        state.normalGap = state.gap(1:3:end);
        state.tangentialSlip = rawSlip;
        state.tangentialGap = stateOld.tangentialGap + rawSlip;
      end

      setState(obj,state);

    end

    function computeContactMatricesAndRhs(obj)

      % A local law returns r, G=dR/dg and Q=dR/dt using the same branch
      % scaling for smooth and semi-smooth contact. Mortar assembly is shared.

      m = MortarSide.master;
      s = MortarSide.slave;
      surfMaster = obj.grids(m).surfaces;
      surfSlave = obj.grids(s).surfaces;
      dofMaster = getDoFManager(obj,m);
      dofSlave = getDoFManager(obj,s);
      fldM = dofMaster.getVariableId(obj.coupledVariables);
      fldS = dofSlave.getVariableId(obj.coupledVariables);
      topolMaster = getRowsMatrix(surfMaster.connectivity,1:surfMaster.num);
      topolSlave = getRowsMatrix(surfSlave.connectivity,1:surfSlave.num);
      [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj);
      rhsUm = zeros(getNumbDoF(dofMaster,obj.coupledVariables),1);
      rhsUs = zeros(getNumbDoF(dofSlave,obj.coupledVariables),1);
      rhsT = zeros(getNumbDoF(obj),1);

      state = getState(obj);
      stateOld = getStateOld(obj);
      stateIni = getStateInit(obj);
      deltaGap = state.gap - stateOld.gap;
      deltaTangentialGap = state.tangentialGap - stateOld.tangentialGap;

      [cN,cT] = getComplementarityParameters(obj);
      residual = zeros(3,surfSlave.num);
      gapTangent = zeros(3,3,surfSlave.num);
      tractionTangent = zeros(3,3,surfSlave.num);
      H = [];

      if obj.isSmooth
        updateAssemblyContactModes(obj,state,deltaTangentialGap);
        for is = 1:surfSlave.num
          id = getMultiplierDoF(obj,is);
          g = [state.gap(3*is-2); deltaGap(3*is-1:3*is)];
          slip = deltaTangentialGap(2*is-1:2*is);
          [residual(:,is),gapTangent(:,:,is),tractionTangent(:,:,is)] = ...
            getLocalContactLaw(obj,obj.activeSet.curr(is),...
            state.traction(id),g,slip,cN,cT);
        end
      else
        % Recompute raw gaps from the current displacement state, so branch
        % selection never depends on the ordering of updateState calls.
        computeGap(obj);
        state = getState(obj);
        deltaGap = state.gap-stateOld.gap;
        rawConstraintGap = deltaGap;
        rawConstraintGap(1:3:end) = state.gap(1:3:end);

        if isempty(obj.stabilizationMat)
          computeStabilizationMatrix(obj);
        end
        % Use a fixed operator: projection derivatives automatically disable
        % stabilization in open directions. No active-set-dependent masking.
        H = obj.stabilizationMat;
        referenceChange = state.traction-stateOld.traction;
        referenceChange(1:3:end) = ...
          state.traction(1:3:end)-stateIni.traction(1:3:end);
        areaSlave = repelem(obj.getSlaveArea(),3,1);
        defect = rawConstraintGap-(H*referenceChange)./areaSlave;

        for is = 1:surfSlave.num
          id = getMultiplierDoF(obj,is);
          [residual(:,is),gapTangent(:,:,is),tractionTangent(:,:,is),mode] = ...
            getSemismoothContactLaw(obj,state.traction(id),defect(id),...
            cN,cT,isForceStickElement(obj,is));
          obj.activeSet.curr(is) = mode;
        end
      end

      elemPairs = obj.quadrature.interfacePairs;
      for vtkSlave = surfSlave.vtkTypes
        elSlave = getElement(obj,vtkSlave,s);
        for vtkMaster = surfMaster.vtkTypes
          elMaster = getElement(obj,vtkMaster,m);
          for iPair = 1:obj.quadrature.numbInterfacePairs
            is = elemPairs(iPair,s);
            im = elemPairs(iPair,m);
            if surfSlave.VTKType(is) ~= vtkSlave; continue; end
            if surfMaster.VTKType(im) ~= vtkMaster; continue; end

            nodesS = surfSlave.loc2glob(topolSlave(is,1:elSlave.nNode));
            nodesM = surfMaster.loc2glob(topolMaster(im,1:elMaster.nNode));
            usDof = dofSlave.getLocalDoF(fldS,nodesS);
            umDof = dofMaster.getLocalDoF(fldM,nodesM);
            tDof = getMultiplierDoF(obj,is);
            [Aum,Aus,area] = getContactPairOperators(obj,iPair,...
              im,is,elMaster,elSlave);
            dTrac = state.traction(tDof) - stateIni.traction(tDof);

            % Equilibrium: displacement-test jump paired with traction.
            asbMu.localAssembly(umDof,tDof,Aum);
            asbDu.localAssembly(usDof,tDof,-Aus);
            rhsUm(umDof) = rhsUm(umDof) + Aum*dTrac;
            rhsUs(usDof) = rhsUs(usDof) - Aus*dTrac;

            % Constraint: one common assembly for stick, slip and open.
            G = gapTangent(:,:,is);
            Q = tractionTangent(:,:,is);
            asbMt.localAssembly(tDof,umDof,G*Aum');
            asbDt.localAssembly(tDof,usDof,-G*Aus');
            asbQ.localAssembly(tDof,tDof,area*Q);
            rhsT(tDof) = rhsT(tDof) + area*residual(:,is);
          end
        end
      end

      obj.addJum(m,asbMu.sparseAssembly());
      obj.addJum(s,asbDu.sparseAssembly());
      obj.addJmu(m,asbMt.sparseAssembly());
      obj.addJmu(s,asbDt.sparseAssembly());
      obj.Jconstraint = asbQ.sparseAssembly();
      if ~obj.isSmooth
        % Chain rule for defect = rawGap - A^{-1} H*(t-reference).
        % A commutes with each face-local 3x3 block; retain ALL off-face
        % traction couplings, including those in sliding rows.
        n = getNumbDoF(obj);
        row = zeros(9*surfSlave.num,1);
        col = row;
        val = row;
        for is = 1:surfSlave.num
          id = getMultiplierDoF(obj,is);
          [rr,cc] = ndgrid(id,id);
          k = (is-1)*9+(1:9);
          row(k) = rr(:);
          col(k) = cc(:);
          block = gapTangent(:,:,is);
          val(k) = block(:);
        end
        Gglobal = sparse(row,col,val,n,n);
        obj.Jconstraint = obj.Jconstraint-Gglobal*H;
      end
      obj.addRhs(m,rhsUm);
      obj.addRhs(s,rhsUs);
      obj.rhsConstraint = rhsT;
    end

    function [Aum,Aus,area] = getContactPairOperators(obj,iPair,...
        im,is,elMaster,elSlave)
      % Pure mortar geometry, independent of contact mode and strategy.
      xiMaster = obj.quadrature.getMasterGPCoords(iPair);
      xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
      dJw = obj.quadrature.getIntegrationWeights(iPair);
      area = sum(dJw);
      [Ns,Nm,Nmult] = getMortarBasisFunctions(obj.quadrature,...
        im,is,elMaster,elSlave,xiMaster,xiSlave);
      [Ns,Nm,Nmult] = reshapeBasisFunctions(3,Ns,Nm,Nmult);
      f = @(a,b) pagemtimes(a,'ctranspose',b,'none');
      R = getRotationMatrix(obj,MortarSide.slave,is);
      Aum = MortarQuadrature.integrate(f,Nm,Nmult,dJw)*R;
      Aus = MortarQuadrature.integrate(f,Ns,Nmult,dJw)*R;
    end

    function updateAssemblyContactModes(obj,~,~)
      % Smooth assembly retains the outer active set. Semi-smooth assembly
      % selects the branch together with its residual and generalized tangent.
      enforceForcedStick(obj);
    end

    function [r,G,Q,mode] = getSemismoothContactLaw(obj,t,d,cN,cT,forced)
      % Projection residual with gap units in ALL branches:
      % R_N = (min(t_N+cN*d_N,0)-t_N)/cN;
      % R_T = (Proj_ball(t_T+cT*d_T)-t_T)/cT.
      % d is an algebraic constraint defect, not the physical gap.
      r = zeros(3,1);
      G = zeros(3);
      Q = zeros(3);
      if forced
        mode = ContactMode.stick;
        r = d;
        G = eye(3);
        return
      end
      qN = t(1)+cN*d(1);
      if qN > 0
        mode = ContactMode.open;
        Q = -diag([1/cN,1/cT,1/cT]);
        r = Q*t;
        return
      end
      r(1) = d(1);
      G(1,1) = 1;
      % Use the projected normal predictor in the Coulomb radius. At a
      % converged closed constraint qN=tN. This also differentiates normal
      % gap and its stabilization coupling in the slipping equation.
      [tauLim,dTauDqN] = getFrictionLimit(obj,qN);
      qT = t(2:3)+cT*d(2:3);
      qNorm = norm(qT);
      if qNorm <= tauLim
        mode = ContactMode.stick;
        r(2:3) = d(2:3);
        G(2:3,2:3) = eye(2);
      else
        mode = ContactMode.slip;
        % qNorm>tauLim>=0: normalization is safe without a cutoff that
        % would invalidate the residual/Jacobian for small sliding loads.
        direction = qT/qNorm;
        Ddirection = (eye(2)-direction*direction')/qNorm;
        r(2:3) = (tauLim*direction-t(2:3))/cT;
        G(2:3,1) = (cN/cT)*dTauDqN*direction;
        G(2:3,2:3) = tauLim*Ddirection;
        Q(2:3,1) = dTauDqN*direction/cT;
        Q(2:3,2:3) = (tauLim*Ddirection-eye(2))/cT;
      end
    end

    function [r,G,Q] = getLocalContactLaw(obj,mode,t,g,slip,cN,cT)
      % Preserve the pasted class branch convention in BOTH strategies:
      % closed normal: gN; stick: dgT; slip: tT-tau*n; open: C\t.
      % The residual and both tangents always use the same row scaling.
      % r: unintegrated residual; G: derivative w.r.t. raw gap/slip;
      % Q: derivative w.r.t. traction. Stabilization is added separately.
      r = zeros(3,1);
      G = zeros(3);
      Q = zeros(3);
      if mode == ContactMode.open
        % Correct the pasted open residual/tangent mismatch: both use 1/c.
        Q = diag([1/cN,1/cT,1/cT]);
        r = Q*t;
        return
      end

      % Closed normal contact is identical in both formulations.
      r(1) = g(1);
      G(1,1) = 1;
      if mode == ContactMode.stick
        r(2:3) = g(2:3);
        G(2:3,2:3) = eye(2);
      elseif mode == ContactMode.slip || mode == ContactMode.newSlip
        % Augmentation enters the Coulomb direction in BOTH strategies.
        [r(2:3),G(2:3,2:3),Q(2:3,:)] = ...
          getCoulombResidualAndTangent(obj,t,slip,cT);
      else
        error('%s: unsupported contact mode.',class(obj));
      end
    end

    function [r,G,Q] = getCoulombResidualAndTangent(obj,t,slip,cT)

      % Unscaled tangential traction residual, as in the pasted class.
      [tauLim,dTauDtN] = getFrictionLimit(obj,t(1));
      [n,Dn] = getUnitVectorAndDerivative(obj,t(2:3) + cT*slip);
      r = t(2:3) - tauLim*n;
      G = -tauLim*cT*Dn;
      Q = [-dTauDtN*n, eye(2)-tauLim*Dn];
      
    end

    function [tauLim,dTauDtN] = getFrictionLimit(obj,tN)
      tanPhi = tan(deg2rad(obj.phi));
      tauRaw = obj.cohesion - tanPhi*tN;
      tauLim = max(tauRaw,0);
      dTauDtN = -tanPhi*double(tauRaw > 0);
    end

    function enforceForcedStick(obj)
      for is = 1:numel(obj.activeSet.curr)
        if isForceStickElement(obj,is)
          obj.activeSet.curr(is) = ContactMode.stick;
        end
      end
    end

    function [cN,cT] = getComplementarityParameters(obj)

      c = obj.contactAugmentation;
      validateattributes(c,{'numeric'},...
        {'vector','numel',2,'real','finite','positive'});

      cN = c(1);
      cT = c(end);

      
    end

    function [n,DnDx] = getUnitVectorAndDerivative(obj,x)
      % Derivative of the normalized trial traction. Retain the existing
      % zero direction/tangent below the sliding tolerance.

      xNorm = norm(x);
      tol = obj.activeSet.tol.sliding;

      if xNorm > tol
        n = x/xNorm;
        DnDx = (eye(numel(x)) - n*n')/xNorm;
      else
        n = zeros(size(x));
        DnDx = zeros(numel(x),numel(x));
      end

    end

    function isForced = isForceStickElement(obj,is)
      % Check whether the current slave surface is constrained to remain in
      % stick mode because of the optional forceStick settings.

      isForced = obj.forceStick;

      if isForced
        return
      end

      if ~isstring(obj.activeSet.forceStickBoundary)
        return
      end

      surfSlave = obj.grids(MortarSide.slave).surfaces;
      nodes = getRowsMatrix(surfSlave.connectivity,is);
      nodes = surfSlave.loc2glob(nodes);

      isForced = any(ismember(nodes,obj.stickNodes));

    end

    function [H,rhsH] = getStabilizationMatrixAndRhs(obj)
      % Keep all stick components and only the normal slip component.
      % Open faces and tangential slip components require no stabilization.
      % H is used directly in the closed/stick gap equations, as in the
      % pasted class. Slip tangential rows and open rows are filtered out.

      if isempty(obj.stabilizationMat)
        computeStabilizationMatrix(obj);
      end

      state = getState(obj);
      iniTrac = getStateInit(obj,"traction");

      H = obj.stabilizationMat;

      elOpen = find(obj.activeSet.curr == ContactMode.open);
      elSlip = [find(obj.activeSet.curr == ContactMode.slip);...
        find(obj.activeSet.curr == ContactMode.newSlip)];

      dofOpen = DoFManager.dofExpand(elOpen,3);
      dofSlip = [3*elSlip-1; 3*elSlip];

      % remove rows and columns of dofs not requiring stabilization
      H([dofOpen;dofSlip],:) = 0;
      H(:,[dofOpen;dofSlip]) = 0;

      % use traction variation for tangential components
      rhsH = -H*state.deltaTraction;

      rhsH(1:3:end) = -H(1:3:end,:) * (state.traction - iniTrac);

    end

    function [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj)
      % helper to define contact matrix assemblers

      s = MortarSide.slave;
      m = MortarSide.master;

      surfMaster = obj.grids(m).surfaces;
      surfSlave = obj.grids(s).surfaces;
      dofSlave = getDoFManager(obj,s);
      dofMaster = getDoFManager(obj,m);

      ncomp = 3;

      elemPairs = obj.quadrature.interfacePairs;
      nv = surfMaster.numVerts(elemPairs(:,m));
      nNMPS = accumarray(elemPairs(:,s),nv,[surfSlave.num,1]);

      N1 = sum(nNMPS);
      N2 = sum(surfSlave.numVerts(elemPairs(:,s)));

      nmu = (ncomp^2)*N1;
      nsu = ncomp^2*N2;
      nmt = nmu;
      nst = nsu;
      nq = ncomp^2*N2;

      nDofMaster = dofMaster.getNumbDoF(obj.coupledVariables);
      nDofSlave = dofSlave.getNumbDoF(obj.coupledVariables);
      nDofMult = getNumbDoF(obj);

      % initialize sparse matrix assemblers
      asbMu = assembler(nmu,nDofMaster,nDofMult);
      asbDu = assembler(nsu,nDofSlave,nDofMult);
      asbMt = assembler(nmt,nDofMult,nDofMaster);
      asbDt = assembler(nst,nDofMult,nDofSlave);
      asbQ = assembler(nq,nDofMult,nDofMult);

    end

    function setStickNodes(obj)

      % set boundary nodes that must remain stick

      bcs = obj.domains(2).bcs;
      bcList = keys(bcs.db);

      if ~isstring(obj.activeSet.forceStickBoundary)
        return
      end

      directions = ismember(["x","y","z"],obj.activeSet.forceStickBoundary);

      stickList = [];

      for bcId = string(bcList)

        if strcmpi(getType(bcs,bcId),"dirichlet") && getVariable(bcs,bcId) == obj.coupledVariables

          nEnts = getNumbTargetEntities(bcs,bcId);

          if sum(nEnts(directions))==0
            continue
          end

          stickList = [stickList; getTargetEntities(bcs,bcId)];

        end

      end

      obj.stickNodes = unique(stickList);
      
    end

  end

  methods (Static)

    function var = getCoupledVariables()
      var = Poromechanics.getField();
    end

  end

end
