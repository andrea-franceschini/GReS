classdef SolidMechanicsContactAugmented < MeshTying

  % solid mechanics solver using piece-wise constant multipliers
  % implments semi-smooth newton rewriting constraint inequalities into complementarity function 

  properties
    phi               % friction angle in deg
    cohesion          % cohesion
    contactHelper
    activeSet
    NLIter = 0
    stickNodes     % boundary nodes where contact state should stay stick
    forceStick     % flag to enforce interface to stay stick
    count
    contactAugmentation
  end


  methods
    function obj = SolidMechanicsContactAugmented(id,domains,inputStruct)

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
        'stabilizationScale',1.0,'augmentationParameter',1.0),input);
      params = readInput(struct('cohesion',[],'frictionAngle',[]),input.Coulomb);

      obj.stabilizationScale = input.stabilizationScale;
      obj.contactAugmentation = input.augmentationParameter;


      obj.forceStick = logical(input.forceStick);

      obj.cohesion = params.cohesion;
      obj.phi = params.frictionAngle;

      nDofsInterface = getNumbDoF(obj);

      s = getState(obj);

      s.traction = zeros(nDofsInterface,1);
      s.deltaTraction = zeros(nDofsInterface,1);

      % the gap in global coordinates
      s.gap = zeros(nDofsInterface,1);
      s.normalGap = zeros(round(1/3*nDofsInterface),1);
      s.tangentialGap = zeros(round(2/3*nDofsInterface),1);
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

      %applyContactReturnMap(obj);

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
      % 

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
      % Project the trial contact multiplier onto the admissible normal cone
      % and Coulomb disk. This implements the return-map stage of the
      % generalized Newton method with return map.
      %
      % Sign convention: normal compression is t_N <= 0. The admissible
      % normal set is therefore t_N <= 0. The tangential admissible set is
      % ||t_T|| <= tau_max, with tau_max = cohesion - tan(phi)*t_N.

      state = getState(obj);
      stateOld = getStateOld(obj);
      surfSlave = obj.grids(MortarSide.slave).surfaces;
      tanPhi = tan(deg2rad(obj.phi));
      %tolTau = obj.activeSet.tol.tangentialViolation;

      for is = 1:numel(obj.activeSet.curr)

        id = DoFManager.dofExpand(is,3);
        tTrial = state.traction(id);

        % Normal return map: t_N = min(t_N^trial,0).
        tN = min(tTrial(1),0.0);

        % If the point is geometrically open and the projected normal
        % traction is zero, release also the tangential traction. This avoids
        % carrying frictional shear on an open interface.
        if state.normalGap(is) > obj.activeSet.tol.normalGap && abs(tN) <= eps
          state.traction(id) = 0;
          obj.activeSet.curr(is) = ContactMode.open;
          continue
        end

        % Tangential return map onto the Coulomb disk.
        tauLim = max(obj.cohesion - tanPhi*tN,0.0);
        tTTrial = tTrial(2:3);
        tTNorm = norm(tTTrial);

        if tTNorm <= tauLim
          tT = tTTrial;
          obj.activeSet.curr(is) = ContactMode.stick;
        else
          if tTNorm > 0
            tT = tauLim*(tTTrial/tTNorm);
            obj.activeSet.curr(is) = ContactMode.slip;
          else
            tT = zeros(2,1);
            obj.activeSet.curr(is) = ContactMode.stick;
          end
        end

        state.traction(id) = [tN; tT];

      end

      state.deltaTraction = state.traction - stateOld.traction;
      setState(obj,state);

    end



    function hasConfigurationChanged = updateConfiguration(obj)

      % Contact status changes are handled inside the semi-smooth Newton
      % residual through complementarity functions. No outer active-set
      % update is required.
      hasConfigurationChanged = false;
      obj.NLIter = 0;

    end


    function initialize(obj)

      initialize@InterfaceSolver(obj);

      % initial traction from cell stress
      tIni = computeInitialTraction(obj);

      addInitialTraction(obj,tIni);

      setStateOld(obj,getState(obj,"traction"),"traction");

      setStickNodes(obj);

    end

    function timeStepSetup(obj)

      % Keep the previous contact classification only as an initial guess for
      % the generalized derivative. The semi-smooth residual updates the
      % diagnostic state during assembly, so no trial-stick reset is applied.
      obj.NLIter = 0;

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


    function computeGap(obj)
      % compute normal gap and tangential slip (local coordinates)

      state = getState(obj);
      stateOld = getStateOld(obj);

      um = obj.domains(MortarSide.master).getState("displacements");
      us = obj.domains(MortarSide.slave).getState("displacements");

      % recover variationally consistent stabilized gaps
      areaSlave = repelem(obj.getSlaveArea(),3,1);

      areaGap = (obj.D*us + obj.M*um);

      state.gap = areaGap./areaSlave;

      [~,rhsStab] = getStabilizationMatrixAndRhs(obj);

      stabGap = (areaGap + rhsStab)./areaSlave;
      stabSlip = (state.gap-stateOld.gap) + rhsStab./areaSlave;

      stabSlip(1:3:end) = [];

      state.tangentialSlip = stabSlip;
      %
      state.normalGap = stabGap(1:3:end);
      state.tangentialGap = stateOld.tangentialGap + stabSlip;

      setState(obj,state);

    end


    function computeContactMatricesAndRhs(obj)

      % Compute contact matrices and rhs using complementarity functions.
      %
      % Normal contact is written in primal-dual form as
      %
      %   min(0, t_N + c_N g_N) - t_N = 0,
      %
      % with the sign convention t_N <= 0 in compression and g_N >= 0 in
      % separation. The assembled residual is scaled by 1/c_N so that the
      % closed-contact branch reduces to the standard mortar constraint
      % g_N = 0.
      %
      % Friction is written in the form used in semi-smooth mortar contact:
      %
      %   ||t_T + c_T dg_T|| <= tau_max  ->  dg_T = 0,
      %   ||t_T + c_T dg_T|| >  tau_max  ->  t_T - tau_max n_T = 0,
      %
      % where n_T = (t_T + c_T dg_T)/||t_T + c_T dg_T|| and
      % tau_max = cohesion - tan(phi) t_N. The sets are therefore not
      % updated by an external configuration loop; they are selected locally
      % by the current Newton iterate and define a generalized derivative.

      prevAS = obj.activeSet.curr;

      m = MortarSide.master;
      s = MortarSide.slave;

      surfMaster = obj.grids(m).surfaces;
      surfSlave = obj.grids(s).surfaces;

      dofMaster = getDoFManager(obj,m);
      dofSlave =  getDoFManager(obj,s);

      elemPairs = obj.quadrature.interfacePairs;

      % define matrix assemblers
      [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj);

      % define rhs vectors
      rhsUm = zeros(getNumbDoF(dofMaster,obj.coupledVariables),1);
      rhsUs = zeros(getNumbDoF(dofSlave,obj.coupledVariables),1);
      rhsT = zeros(getNumbDoF(obj),1);

      fldM = dofMaster.getVariableId(obj.coupledVariables);
      fldS = dofSlave.getVariableId(obj.coupledVariables);

      % anonymous functions for local mortar computations
      f1 = @(a,b) pagemtimes(a,'ctranspose',b,'none');

      state = getState(obj);
      stateIni = getStateInit(obj);

      % Use the incremental tangential gap as the frictional slip variable.
      stateOld = getStateOld(obj);
      deltaGap = state.gap - stateOld.gap;
      deltaTangentialGap = state.tangentialGap - stateOld.tangentialGap;

      topolMaster = getRowsMatrix(surfMaster.connectivity,1:surfMaster.num);
      topolSlave = getRowsMatrix(surfSlave.connectivity,1:surfSlave.num);

      tols = obj.activeSet.tol;

      % The state is now diagnostic. It is still useful for VTK output,
      % stabilization filtering and debugging.
      obj.activeSet.curr(:) = ContactMode.stick;

      for vtkSlave = surfSlave.vtkTypes

        elSlave = getElement(obj,vtkSlave,s);

        for vtkMaster = surfMaster.vtkTypes

          elMaster = getElement(obj,vtkMaster,m);

          % loop over pairs of connected master/slave elements
          for iPair = 1:obj.quadrature.numbInterfacePairs

            is = elemPairs(iPair,s);
            im = elemPairs(iPair,m);

            if surfSlave.VTKType(is) ~= vtkSlave; continue; end
            if surfMaster.VTKType(im) ~= vtkMaster; continue; end

            % retrieve mortar integration data
            xiMaster = obj.quadrature.getMasterGPCoords(iPair);
            xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
            dJw = obj.quadrature.getIntegrationWeights(iPair);

            % area of current integration cell
            area = sum(dJw);

            % define slave related quantities
            nodeSlave = surfSlave.loc2glob(topolSlave(is,1:elSlave.nNode));
            usDof = dofSlave.getLocalDoF(fldS,nodeSlave);
            tDof = getMultiplierDoF(obj,is);
            trac = state.traction(tDof);
            tIni = stateIni.traction(tDof);

            % equilibrium equation and stabilization work with the traction
            % variation, so that an initial balanced traction does not add
            % spurious internal work.
            dTrac = trac - tIni;

            nodeMaster = surfMaster.loc2glob(topolMaster(im,1:elMaster.nNode));
            umDof = dofMaster.getLocalDoF(fldM,nodeMaster);

            [Nslave,Nmaster,Nmult] = ...
              getMortarBasisFunctions(obj.quadrature,im,is,elMaster,elSlave,xiMaster,xiSlave);

            % reshape basis function matrices to match number of components
            [Ns,Nm,Nmult] = reshapeBasisFunctions(3,Nslave,Nmaster,Nmult);

            % rotation matrix
            R = getRotationMatrix(obj,MortarSide.slave,is);

            % local gap variables
            g_n = state.gap(3*is-2);
            dgt = deltaGap([3*is-1; 3*is]);
            dgtStab = deltaTangentialGap([2*is-1; 2*is]);

            % A_us
            Aum =  MortarQuadrature.integrate(f1,Nm,Nmult,dJw);
            Aus =  MortarQuadrature.integrate(f1,Ns,Nmult,dJw);

            % apply rotation matrix due to mixed dof assembly
            Aum = Aum*R;
            Aus = Aus*R;

            asbMu.localAssembly(umDof,tDof,Aum);
            asbDu.localAssembly(usDof,tDof,-Aus);

            % rhs (jump(eta),t)
            rhsUm(umDof) = rhsUm(umDof) + Aum*dTrac;
            rhsUs(usDof) = rhsUs(usDof) - Aus*dTrac;

            % Derivatives of the local gap with respect to displacement
            BgN_m = Aum(:,1)';
            BgN_s = Aus(:,1)';
            BgT_m = Aum(:,2:3)';
            BgT_s = Aus(:,2:3)';

            [cN,cT] = getComplementarityParameters(obj,area);
            tanPhi = tan(deg2rad(obj.phi));

            tN = trac(1);
            tT = trac(2:3);

            zN = tN + cN*g_n;
            tauLim = max(obj.cohesion - tanPhi*tN,0.0);


            % --- normal complementarity ---------------------------------
            if zN > 0

              % open branch: t_N = 0 and t_T = 0
              obj.activeSet.curr(is) = ContactMode.open;

              rhsT(tDof(1)) = rhsT(tDof(1)) + area*tN/cN;
              asbQ.localAssembly(tDof(1),tDof(1),area/cN);

              rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + area*tT/cT;
              asbQ.localAssembly(tDof(2:3),tDof(2:3),area/cT*eye(2));

              continue

            else

              % closed branch: g_N = 0
              rhsT(tDof(1)) = rhsT(tDof(1)) + area*g_n;
              asbMt.localAssembly(tDof(1),umDof,BgN_m);
              asbDt.localAssembly(tDof(1),usDof,-BgN_s);

            end

            contactState = obj.activeSet.curr(is);

            % --- tangential complementarity ------------------------------
            cT = min([norm(tauLim)/norm(dgtStab),1e4]);
            yT = tT + cT*dgtStab;     % trial tangential traction
            yNorm = norm(yT);
            tau = yNorm;

            % apply hysteresis to tangential traction norm for check
            if contactState == ContactMode.stick && tau >= tauLim

              % reduce the tau if goes above limit
              tau = tau*(1-tols.tangentialViolation);

            elseif contactState ~= ContactMode.stick  && tau <=tauLim

              % increase tau if falls below limit
              tau = tau*(1+tols.tangentialViolation);
            end


            if tau <= tauLim

              % stick branch: dg_T = 0
              obj.activeSet.curr(is) = ContactMode.stick;

              rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + area*dgt;
              asbMt.localAssembly(tDof(2:3),umDof,BgT_m);
              asbDt.localAssembly(tDof(2:3),usDof,-BgT_s);

            else

              % slip branch: t_T = tau_max * n_T
              obj.activeSet.curr(is) = ContactMode.slip;

              [nT,DnDy] = getUnitVectorAndDerivative(obj,yT);

              RT = tT - tauLim*nT;

              rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + area*RT/cT;

              % dR_T/dt_N = -d(tauLim)/d(tN) * n_T
              tauRaw = obj.cohesion - tanPhi*tN;
              if tauRaw > 0
                dTauDtN = -tanPhi;
              else
                dTauDtN = 0;
              end
              asbQ.localAssembly(tDof(2:3),tDof(1),-area*dTauDtN*nT/cT);


              % dR_T/dt_T = I - tau_max d(n_T)/d(y_T)
              dRdtT = eye(2) - tauLim*DnDy;
              asbQ.localAssembly(tDof(2:3),tDof(2:3),area*dRdtT/cT);

              % dR_T/ddg_T = -tau_max d(n_T)/d(y_T) c_T
              dRdgt = -tauLim*DnDy;
              asbMt.localAssembly(tDof(2:3),umDof,dRdgt*BgT_m);
              asbDt.localAssembly(tDof(2:3),usDof,-dRdgt*BgT_s);

            end

          end % end inner master elems loop

        end
      end

      % if ~any(int32(obj.activeSet.curr) - int32(prevAS))
      %   fprintf('Active set did not change \n')
      % end


      % assemble matrices into jacobian blocks
      obj.addJum(MortarSide.master, asbMu.sparseAssembly());
      obj.addJum(MortarSide.slave, asbDu.sparseAssembly());
      obj.addJmu(MortarSide.master, asbMt.sparseAssembly());
      obj.addJmu(MortarSide.slave, asbDt.sparseAssembly());

      obj.Jconstraint = asbQ.sparseAssembly();

      obj.addRhs(MortarSide.master,rhsUm);
      obj.addRhs(MortarSide.slave,rhsUs);
      obj.rhsConstraint = rhsT;

    end


    function [cN,cT] = getComplementarityParameters(obj,area)
      % Return the local augmentation parameters used by the primal-dual
      % complementarity functions. A scalar input is used for both normal
      % and tangential directions. A two-entry input can be used to prescribe
      % different normal and tangential values.

      c = obj.contactAugmentation;

      if isempty(c)
        c = 1.0;
      end

      if isscalar(c)
        cN = c;
        cT = c;
      else
        cN = c(1);
        cT = c(2);
      end

      if cN <= 0 || cT <= 0
        error('%s: augmentationParameter must be strictly positive.',class(obj));
      end

    end


    function [n,DnDx] = getUnitVectorAndDerivative(obj,x)
      % Generalized derivative of x/||x||. At the origin, pick a bounded
      % element of the generalized derivative to avoid division by zero.

      xNorm = norm(x);
      tol = obj.activeSet.tol.sliding;

      if xNorm > tol
        n = x/xNorm;
        DnDx = (eye(numel(x)) - n*n')/xNorm;
      else
        n = zeros(size(x));
        DnDx = zeros(numel(x));
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
      % reutnr the stabilization matrix after removing contribution for traction
      % dofs that do not need stabilization:
      % - normal component in slip dofs
      % - all components of open dofs

      if isempty(obj.stabilizationMat)
        computeStabilizationMatrix(obj);
      end

      state = getState(obj);
      iniTrac = getStateInit(obj,"traction");

      H = obj.stabilizationMat;

      %
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

    function dtdgt = computeDerTracGap(obj,sigma_n,slip)
      % gt = obj.g_T(get_dof(nodeId));

      % slip: 2x1 local tangential slip
      tauLim = obj.cohesion - tan(deg2rad(obj.phi))*sigma_n;
      dtdgt = tauLim*((eye(2)*norm(slip)^2 - slip*slip')/(norm(slip))^3);

    end

    function dtdtn = computeDerTracTn(obj,slip,t)

      tanPhi = tan(deg2rad(obj.phi));

      if norm(slip) > obj.activeSet.tol.sliding
        %use available gap to properly compute traction
        dtdtn = -tanPhi*(slip/norm(slip));
      else
        t = t(2:3);
        dtdtn = -tanPhi*(t/norm(t));
      end
    end

    function tracLim = computeLimitTraction(obj,dgt,t,slipNorm)

      % return the limit traction vector in the local frame
      t_N = t(1);
      tauLim = obj.cohesion - tan(deg2rad(obj.phi))*t_N;


      if slipNorm > obj.activeSet.tol.sliding && obj.NLIter > 0
        tracLim = tauLim*(dgt/norm(dgt));
      else
        % compute tangential traction from traction (global coordinates!)
        t = t(2:3);
        tracLim =  tauLim*(t/norm(t));
      end

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


