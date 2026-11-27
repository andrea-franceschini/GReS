classdef SolidMechanicsContact < MortarContinuity

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
    forcedStickCells
    contactAugmentation
    isSmooth = true % outer active-set loop; false selects semi-smooth Newton
  end

  methods

    function obj = SolidMechanicsContact(interface)

      obj@MortarContinuity(interface);

    end

    function registerConstraint(obj,varargin)

      registerConstraint@MortarContinuity(obj,varargin{:});

      if ~strcmp(obj.multiplierType,'P0')
        error("Constraint %s is implemented only for P0 multipliers" + ...
          " located at %s",...
          class(obj),entityField.cell)
      end


      input = readInput(struct('Coulomb',[],'ActiveSet',missing,'forceStick',0,'augmentationParameter',1.0,...
        'augmentationNormal',missing,'augmentationTangential',missing,'isSmooth',true),varargin{:});

      params = readInput(struct('cohesion',0.0,'frictionAngle',[]),input.Coulomb);


      if ~ismissing(input.augmentationNormal)
        obj.contactAugmentation(1) = input.augmentationNormal;
      else
        obj.contactAugmentation(1) = input.augmentationParameter;
      end

      if ~ismissing(input.augmentationTangential)
        obj.contactAugmentation(2) = input.augmentationTangential;
      else
        obj.contactAugmentation(2) = input.augmentationParameter;
      end

      obj.isSmooth = logical(input.isSmooth);

      obj.forceStick = logical(input.forceStick);

      obj.cohesion = params.cohesion;
      obj.phi = params.frictionAngle;

      dofm = obj.interface.dofm;

      fld = obj.getField;

      nTracDof = dofm.getNumbDoF(fld);

      s = getState(obj);

      s.traction = zeros(nTracDof,1);
      s.deltaTraction = zeros(nTracDof,1);

      % Raw mortar gap in the local contact frame
      s.gap = zeros(nTracDof,1);
      s.normalGap = zeros(round(1/3*nTracDof),1);

      % total tangential gap
      s.tangentialGap = zeros(round(2/3*nTracDof),1);

      % gap variation at the current time step
      s.tangentialSlip = zeros(round(2/3*nTracDof),1);
      setState(obj,s);

      initializeActiveSet(obj,obj.grid.cells.num,input.ActiveSet);

    end

    function updateState(obj,du)

      % traction update
      actMult = getMultiplierDoF(obj);

      state = getState(obj);
      stateOld = getStateOld(obj);

      dofm = obj.interface.dofm;
      dofs = dofm.getDoF(dofm.getVariableId(obj.getField));
      state.traction(actMult) = state.traction(actMult) + du(dofs);
      
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

      interf = obj.interface;
      model = interf.model;
      master = interf.domains(MortarSide.master);
      slave = interf.domains(MortarSide.slave);
      dispFld = obj.getCoupledField();
      tracFld = obj.getField();

      % clear previous jacobian coupling
      % obj.model.setJacobian(slave,[],dispFld,tracFld,interf);
      % obj.model.setJacobian(master,[],dispFld,tracFld,interf);
      % obj.model.setJacobian(interf,[],tracFld,dispFld,master);
      % obj.model.setJacobian(interf,[],tracFld,dispFld,slave);

      obj.computeConstraintMatrices(dispFld);


      computeContactMatricesAndRhs(obj);

      % get stabilization matrix depending on the current active set
      [H,rhsStab] = getStabilizationMatrixAndRhs(obj);

      model.addJacobian(interf,-H,tracFld);
      model.addResidual(interf,rhsStab,tracFld);

      if gresLog().getVerbosity > 3
        % print rhs terms for each fracture state for debug purposes
        res = model.getResidual(interf,tracFld);
        dof_stick = DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.stick),3);
        dof_slip = [DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.slip),3); ...
          DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.newSlip),3)];
        dof_open = DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.open),3);
        fprintf('Rhs norm for stabilization: %4.3e \n', norm(rhsStab));
        fprintf('Rhs norm for stick dofs: %4.3e \n', norm(res(dof_stick)))
        fprintf('Rhs norm for slip dofs: %4.3e \n', norm(res(dof_slip)))
        fprintf('Rhs norm for open dofs: %4.3e \n', norm(res(dof_open)))
      end

    end


    % function assembleConstraintOld(obj)
    %
    %   % reset the jacobian blocks
    %
    %   interf = obj.interface;
    %   model = interf.model;
    %   master = interf.domains(MortarSide.master);
    %   slave = interf.domains(MortarSide.slave);
    %   dispFld = obj.getCoupledField();
    %   tracFld = obj.getField();
    %
    %   % clear previous jacobian coupling
    %   obj.setJacobian(slave,[],dispFld,tracFld,interf);
    %   obj.setJacobian(master,[],dispFld,tracFld,interf);
    %   obj.setJacobian(interf,[],tracFld,dispFld,master);
    %   obj.setJacobian(interf,[],tracFld,dispFld,slave);
    %
    %   if isempty(obj.Dmat)
    %     computeConstraintMatrices(obj,dispFld);
    %   end
    %
    %   computeContactMatricesAndRhs(obj);
    %
    %   % get stabilization matrix depending on the current active set
    %   [H,rhsStab] = getStabilizationMatrixAndRhs(obj);
    %
    %   model.addJacobian(interf,-H,tracFld);
    %   model.addResidual(interf,rhsStab,tracFld);
    %
    %   if gresLog().getVerbosity > 3
    %     % print rhs terms for each fracture state for debug purposes
    %     res = model.getResidual(interf,tracFld);
    %     dof_stick = DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.stick),3);
    %     dof_slip = [DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.slip),3); ...
    %                 DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.newSlip),3)];
    %     dof_open = DoFManager.dofExpand(find(obj.activeSet.curr == ContactMode.open),3);
    %     fprintf('Rhs norm for stabilization: %4.3e \n', norm(rhsStab));
    %     fprintf('Rhs norm for stick dofs: %4.3e \n', norm(res(dof_stick)))
    %     fprintf('Rhs norm for slip dofs: %4.3e \n', norm(res(dof_slip)))
    %     fprintf('Rhs norm for open dofs: %4.3e \n', norm(res(dof_open)))
    %   end
    %
    % end

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

      initialize@MortarContinuity(obj);

      % initial traction from cell stress
      tIni = computeInitialTraction(obj);

      addInitialTraction(obj,tIni);

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

    end

    function trac = computeInitialTraction(obj)
      quad = obj.interface.quadrature;
      mesh = quad.grid;
      slave = obj.interface.domains(MortarSide.slave);
      slaveFaces = mesh.traceGrids(MortarSide.slave).cells.faceId;
      bulkCells = slave.grid.faces.neighbors(slaveFaces,1);
      avgStress = slave.getState("avgStress");

      % Project the slave stress onto the interface normal, then average
      % traction over the part of the interface cell covered by triplets.
      idx = [1;6;5;6;2;4;5;4;3];
      traction = zeros(3,mesh.cells.num);
      area = zeros(mesh.cells.num,1);
      for triplet = reshape(quad.getActiveTriplets(),1,[])
        cellId = mesh.neighbors(triplet,1);
        slaveId = mesh.neighbors(triplet,2);
        sigma = reshape(avgStress(bulkCells(slaveId),idx),3,3);
        weight = sum(quad.getIntegrationWeights(triplet));
        traction(:,cellId) = traction(:,cellId) + weight*sigma*mesh.cells.normal(cellId,:)';
        area(cellId) = area(cellId) + weight;
      end

      trac = zeros(getNumbDoF(obj),1);
      for cellId = 1:mesh.cells.num
        tDof = getMultiplierDoF(obj,cellId);
        R = obj.getMultiplierRotation(cellId);
        trac(tDof) = R'*traction(:,cellId)/area(cellId);
      end
    end

    function advanceState(obj)

      advanceState@MortarContinuity(obj);
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
      goBackState@MortarContinuity(obj);
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

      state = obj.interface.state;

      outTraction = state.interpolate(fac,"traction");
      dT = state.interpolate(fac,"deltaTraction");
      outNormalGap = state.interpolate(fac,"normalGap");
      outTangentialSlip = state.interpolate(fac,"tangentialSlip");
      outTangentialGap = state.interpolate(fac,"tangentialGap");

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
        'rotationMatrix',          obj.grid.cells.rotationMatrices
        'deltaTraction',           deltaTrac
        };

      surfaceStr = cell2struct(entries, {'name','data'}, 2);
    end

    function writeSolution(obj,fac,tID)

      s = obj.interface.state.interpolate(fac);

      tT = [s.traction(2:3:end),s.traction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      out = obj.interface.outstate;
      out.results(tID).tractionVec = s.traction;
      out.results(tID).normalGap = s.normalGap;
      out.results(tID).slipIncrement = s.tangentialSlip;
      out.results(tID).tangentialGap = s.tangentialGap;
      out.results(tID).tangentialTractionNorm = norm_tT;

    end


    function coupledFld = getCoupledField(~,varargin)

      coupledFld = "displacements";
    end

    function fld = getField(~)

      fld = "traction";
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
      cells = obj.grid.cells;

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

        areaChanged = sum(cells.area(hasChangedElem));
        totArea = sum(cells.area);
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

      um = obj.interface.domains(MortarSide.master).getState("displacements");
      us = obj.interface.domains(MortarSide.slave).getState("displacements");

      % Recover the stabilized gap using the pasted class convention:
      % closed/stick constraints are A*g + rhsStab = 0 in both strategies.
      quad = obj.interface.quadrature;
      % The traction rows are indexed by interface cell, so sum triplet
      % weights on those cells before normalizing the integrated gap.
      area = accumarray(quad.grid.neighbors(:,1),cellfun(@sum,quad.weights), ...
        [quad.grid.cells.num,1]);
      areaDof = zeros(obj.getNumbDoF(),1);
      dofm = obj.getDoFManager();
      fieldId = dofm.getVariableId(obj.getField());
      dofs = dofm.getLocalDoF(fieldId,(1:numel(area)).');
      areaDof(dofs) = repelem(area,3,1);

      [D,M] = obj.computeConstraintMatrices(obj.getCoupledField());

      areaGap = (D*us + M*um);

      state.gap = areaGap./areaDof;

      [~,rhsStab] = getStabilizationMatrixAndRhs(obj);

      stabGap = (areaGap + rhsStab)./areaDof;
      stabSlip = (state.gap-stateOld.gap) + rhsStab./areaDof;

      stabSlip(1:3:end) = [];

      state.tangentialSlip = stabSlip;
      state.normalGap = stabGap(1:3:end);
      state.tangentialGap = stateOld.tangentialGap + stabSlip;

      setState(obj,state);

    end

    function computeContactMatricesAndRhsOld(obj)
      computeContactMatricesAndRhs(obj);
    end

    function computeContactMatricesAndRhs(obj)

      interf = obj.interface;
      model = interf.model;
      m = MortarSide.master;
      s = MortarSide.slave;
      master = interf.domains(m);
      slave = interf.domains(s);
      quad = interf.quadrature;
      mesh = quad.grid;
      meshMaster = mesh.traceGrids(m);
      meshSlave = mesh.traceGrids(s);
      coupFld = obj.getCoupledField();
      tracFld = obj.getField();
      dofMaster = master.dofm;
      dofSlave = slave.dofm;
      fldM = dofMaster.getVariableId(coupFld);
      fldS = dofSlave.getVariableId(coupFld);
      interfaceElements = obj.createTraceElements(mesh,quad.gaussOrder);
      masterElements = obj.createTraceElements(meshMaster,quad.gaussOrder);
      slaveElements = obj.createTraceElements(meshSlave,quad.gaussOrder);

      [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj);
      rhsUm = zeros(dofMaster.getNumbDoF(fldM),1);
      rhsUs = zeros(dofSlave.getNumbDoF(fldS),1);
      rhsT = zeros(getNumbDoF(obj),1);
      state = getState(obj);
      stateOld = getStateOld(obj);
      stateIni = getStateInit(obj);
      deltaGap = state.gap - stateOld.gap;
      deltaTangentialGap = state.tangentialGap - stateOld.tangentialGap;

      % update contact state if semismooth
      updateAssemblyContactModes(obj,state,deltaTangentialGap);
      [cN,cT] = getComplementarityParameters(obj);

      % The contact law belongs to the interface cell. Several trace
      % triplets can contribute to the same P0 traction multiplier.

      % cache tangent matrices
      residual = zeros(3,mesh.cells.num);
      gapTangent = zeros(3,3,mesh.cells.num);
      tractionTangent = zeros(3,3,mesh.cells.num);
      for cellId = 1:mesh.cells.num
        tDof = getMultiplierDoF(obj,cellId);
        g = [state.gap(tDof(1));deltaGap(tDof(2:3))];
        slip = deltaTangentialGap(2*cellId-1:2*cellId);
        [residual(:,cellId),gapTangent(:,:,cellId),tractionTangent(:,:,cellId)] = ...
          getLocalContactLaw(obj,obj.activeSet.curr(cellId),state.traction(tDof),g,slip,cN,cT);
      end

      % loop over facet neighbors
      for triplet = reshape(quad.getActiveTriplets(),1,[])
        ids = mesh.neighbors(triplet,:);
        cellId = ids(1);
        is = ids(s+1);
        im = ids(m+1);
        % get FEM objects
        elemI = interfaceElements{mesh.cells.VTKType(cellId)};
        elemS = slaveElements{meshSlave.cells.VTKType(is)};
        elemM = masterElements{meshMaster.cells.VTKType(im)};

        % get dofs
        nodesS = meshSlave.getCellNodes(is);
        nodesM = meshMaster.getCellNodes(im);
        usDof = dofSlave.getLocalDoF(fldS,meshSlave.cells.loc2glob(nodesS(1:elemS.nNode)));
        umDof = dofMaster.getLocalDoF(fldM,meshMaster.cells.loc2glob(nodesM(1:elemM.nNode)));
        tDof = getMultiplierDoF(obj,cellId);

        % get cross grid contribution for the curretn triplet
        [Aum,Aus,area] = getContactPairOperators(obj,triplet,cellId,elemI,elemM,elemS);

        dTrac = state.traction(tDof) - stateIni.traction(tDof);
        asbMu.localAssembly(umDof,tDof,Aum);
        asbDu.localAssembly(usDof,tDof,-Aus);
        rhsUm(umDof) = rhsUm(umDof) + Aum*dTrac;
        rhsUs(usDof) = rhsUs(usDof) - Aus*dTrac;

        % assemble local tangent blocks
        G = gapTangent(:,:,cellId);
        Q = tractionTangent(:,:,cellId);
        asbMt.localAssembly(tDof,umDof,G*Aum');
        asbDt.localAssembly(tDof,usDof,-G*Aus');
        asbQ.localAssembly(tDof,tDof,area*Q);
        rhsT(tDof) = rhsT(tDof) + area*residual(:,cellId);
      end

      % add jacobian and rhs 
      model.addJacobian(master,asbMu.sparseAssembly(),coupFld,tracFld,interf);
      model.addJacobian(slave,asbDu.sparseAssembly(),coupFld,tracFld,interf);
      model.addJacobian(interf,asbMt.sparseAssembly(),tracFld,coupFld,master);
      model.addJacobian(interf,asbDt.sparseAssembly(),tracFld,coupFld,slave);
      model.addJacobian(interf,asbQ.sparseAssembly(),tracFld);
      model.addResidual(master,rhsUm,coupFld);
      model.addResidual(slave,rhsUs,coupFld);
      model.addResidual(interf,rhsT,tracFld);
    end

    function [Aum,Aus,area] = getContactPairOperators(obj,triplet,cellId,elemI,elemM,elemS)
      quad = obj.interface.quadrature;
      dJw = quad.getIntegrationWeights(triplet);
      area = sum(dJw);
      [Ns,Nm,Nmult] = quad.getMortarBasisFunctions(triplet,elemI,elemM,elemS,obj.multiplierType);
      [Ns,Nm,Nmult] = reshapeBasisFunctions(3,Ns,Nm,Nmult);
      f = @(a,b) pagemtimes(a,'ctranspose',b,'none');
      R = obj.getMultiplierRotation(cellId);
      Aum = MortarQuadrature.integrate(f,Nm,Nmult,dJw)*R;
      Aus = MortarQuadrature.integrate(f,Ns,Nmult,dJw)*R;
    end

    function updateAssemblyContactModes(obj,state,slip)

      if obj.isSmooth
        enforceForcedStick(obj);
        return
      end

      [cN,cT] = getComplementarityParameters(obj);

      for is = 1:numel(obj.activeSet.curr)

        if isForceStickElement(obj,is)
          obj.activeSet.curr(is) = ContactMode.stick;
          continue
        end

        id = getMultiplierDoF(obj,is);
        t = state.traction(id);
        if t(1) + cN*state.normalGap(is) > 0
          obj.activeSet.curr(is) = ContactMode.open;
        else
          tauLim = getFrictionLimit(obj,t(1));
          tau = norm(t(2:3) + cT*slip(2*is-1:2*is));
          % Preserve the tangential violation allowance of Augmented.
          % Its assembly reset each face to stick before this comparison.
          if tau >= tauLim
            tau = tau*(1-obj.activeSet.tol.tangentialViolation);
          end

          if tau <= tauLim
            obj.activeSet.curr(is) = ContactMode.stick;
          else
            obj.activeSet.curr(is) = ContactMode.slip;
          end
        end
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
      isForced = obj.forceStick || ...
        (~isempty(obj.forcedStickCells) && obj.forcedStickCells(is));
    end

    function [H,rhsH] = getStabilizationMatrixAndRhs(obj)
      % Keep all stick components and only the normal slip component.
      % Open faces and tangential slip components require no stabilization.
      % H is used directly in the closed/stick gap equations, as in the
      % pasted class. Slip tangential rows and open rows are filtered out.

      if isempty(obj.stabilizationMat)
        computeStabilizationMatrix(obj,obj.getCoupledField(),obj.getField());
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

      meshMaster = obj.grid.traceGrids(m).cells;
      meshSlave = obj.grid.traceGrids(s).cells;
      dofSlave = getDoFManager(obj,s);
      dofMaster = getDoFManager(obj,m);

      ncomp = 3;

      neigh = obj.grid.neighbors;
      nmu = ncomp^2*sum(meshMaster.numVerts(neigh(:,m+1)));
      nsu = ncomp^2*sum(meshSlave.numVerts(neigh(:,s+1)));
      nmt = nmu;
      nst = nsu;
      nq = ncomp^2*size(neigh,1);

      nDofMaster = dofMaster.getNumbDoF(obj.coupledFields);
      nDofSlave = dofSlave.getNumbDoF(obj.coupledFields);
      nDofMult = getNumbDoF(obj);

      % initialize sparse matrix assemblers
      asbMu = assembler(nmu,nDofMaster,nDofMult);
      asbDu = assembler(nsu,nDofSlave,nDofMult);
      asbMt = assembler(nmt,nDofMult,nDofMaster);
      asbDt = assembler(nst,nDofMult,nDofSlave);
      asbQ = assembler(nq,nDofMult,nDofMult);

    end

    function setStickNodes(obj)
      obj.forcedStickCells = false(obj.grid.cells.num,1);
      obj.stickNodes = [];
      selection = obj.activeSet.forceStickBoundary;
      master = obj.interface.domains(MortarSide.master);
      bcs = master.bcs;
      if isempty(bcs) || ~(isstring(selection) || ischar(selection))
        return
      end

      directions = ismember(["x","y","z"],string(selection));
      stickList = [];
      for bcId = reshape(bcs.getBCList(),1,[])
        if bcs.getType(bcId) == BCtype.dirichlet && ...
            bcs.getVariable(bcId) == obj.getCoupledField()
          counts = bcs.getNumbTargetEntities(bcId);
          nodes = bcs.getTargetEntities(bcId);
          % BC targets are concatenated by component.
          offsets = [0;cumsum(counts(:))];
          for component = find(directions)
            rows = offsets(component)+1:offsets(component+1);
            stickList = [stickList;reshape(nodes(rows),[],1)];
          end
        end
      end
      obj.stickNodes = unique(stickList);

      % Match the selected bulk nodes on the master side, then use the
      % triplets to flag the interface cells supported by those trace facets.
      trace = obj.grid.traceGrids(MortarSide.master);
      isStickNode = ismember(trace.cells.loc2glob,obj.stickNodes);
      isStickFacet = false(trace.cells.num,1);
      for facet = 1:trace.cells.num
        nodes = trace.getCellNodes(facet);
        isStickFacet(facet) = any(isStickNode(nodes(nodes>0)));
      end
      neighbors = obj.grid.neighbors;
      obj.forcedStickCells(neighbors(isStickFacet(neighbors(:,3)),1)) = true;
    end

  end

  methods (Static)

    function var = getCoupledVariables()
      var = Poromechanics.getField();
    end

  end

end
