classdef SolidMechanicsContactRegularized < MeshTying

  % solid mechanics solver using piece-wise constant multipliers

  properties
    phi               % friction angle in deg
    cohesion          % cohesion
    contactHelper
    activeSet
    NLIter = 0
    stickNodes     % boundary nodes where contact state should stay stick
    forceStick     % flag to enforce interface to stay stick
    alphaMin      % minimum weight of traction direction in regularized slip update
    alphaMax
  end


  methods
    function obj = SolidMechanicsContactRegularized(id,domains,inputStruct)

      obj@MeshTying(id,domains,inputStruct);

      if obj.multiplierLocation ~= entityField.surface
        error("Interface Solver %s is not implemented for multipliers" + ...
          "located at %s. The only available entityField is %s",...
          class(obj),obj.multiplierLocation,entityField.surface)
      end

    end

    function registerInterface(obj,varargin)

      input = varargin{1};

      input = readInput(struct('Coulomb',[],'ActiveSet',missing,'forceStick',0,'stabilizationScale',1.0,'alphaMin',0.0,'alphaMax',1.0),input);
      params = readInput(struct('cohesion',[],'frictionAngle',[]),input.Coulomb);

      obj.stabilizationScale = input.stabilizationScale;


      obj.forceStick = logical(input.forceStick);

      obj.cohesion = params.cohesion;
      obj.phi = params.frictionAngle;
      obj.alphaMin = input.alphaMin;
      obj.alphaMax = input.alphaMax;

      % if obj.alphaMin < 0 || obj.alphaMin > 1
      %   error('SolidMechanicsContactRegularized: alphaMin must be in [0,1].')
      % end

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

      % update gap
      computeGap(obj);
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


    function hasConfigurationChanged = updateConfiguration(obj)

      hasConfigurationChanged = false;
      if obj.forceStick
        return
      end

      obj.NLIter = 0;

      oldActiveSet = obj.activeSet.curr;
      surfSlave = obj.grids(MortarSide.slave).surfaces;

      state = getState(obj);

      for is = 1:numel(obj.activeSet.curr)

        currAS = obj.activeSet.curr(is);

        nodes = getRowsMatrix(surfSlave.connectivity,is);
        nodes = surfSlave.loc2glob(nodes);

        if isstring(obj.activeSet.forceStickBoundary)
          % force elements adjacent to dirichlet boundary to remain stick
          if any(ismember(nodes,obj.stickNodes))
            % element has a dirichlet node - keep it stick
            continue
          end
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

      % do not upate state of element that exceeded the maximum number of
      % individual updates
      reset = obj.activeSet.stateChange >= ...
        obj.activeSet.tol.maxStateChange;

      asNew(reset) = asOld(reset);

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
          %obj.activeSet.curr = oldActiveSet;
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


    function initialize(obj)

      initialize@InterfaceSolver(obj);

      % initial traction from cell stress
      tIni = computeInitialTraction(obj);

      addInitialTraction(obj,tIni);

      setStateOld(obj,getState(obj,"traction"),"traction");

      setStickNodes(obj);

    end

    function timeStepSetup(obj)

      % set all active dofs to stick state (force trial stick state)
      isActive = obj.activeSet.curr ~= ContactMode.open;
      %obj.activeSet.curr(isActive) = ContactMode.stick;

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

      % Compute contact matrices and rhs

      % the syntax of the comments is inspired by the appendix of
      % Franceschini, Andrea, et al. "Algebraically stabilized Lagrange
      % multiplier method for frictional contact mechanics with
      % hydraulically active fractures." Computer Methods in Applied
      % Mechanics and Engineering 368 (2020): 113161.

      m = MortarSide.master;
      s = MortarSide.slave;

      surfMaster = obj.grids(m).surfaces;
      surfSlave = obj.grids(s).surfaces;

      dofMaster = getDoFManager(obj,m);
      dofSlave =  getDoFManager(obj,s);

      elemPairs = obj.quadrature.interfacePairs;

      % define matrix assemblers
      [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj);

      % define rhs vectors;
      rhsUm = zeros(getNumbDoF(dofMaster,obj.coupledVariables),1);
      rhsUs = zeros(getNumbDoF(dofSlave,obj.coupledVariables),1);
      rhsT = zeros(getNumbDoF(obj),1);

      fldM = dofMaster.getVariableId(obj.coupledVariables);
      fldS = dofSlave.getVariableId(obj.coupledVariables);

      % anonymous functions for local fem computations
      f1 = @(a,b) pagemtimes(a,'ctranspose',b,'none');
      f2 = @(a,b) pagemtimes(a,b);

      % compute slip
      state = getState(obj);
      stateOld = getStateOld(obj);
      stateIni = getStateInit(obj);
      slip = state.gap - stateOld.gap;
      tangSlip = state.tangentialGap - stateOld.tangentialGap;

      topolMaster = getRowsMatrix(surfMaster.connectivity,1:surfMaster.num);
      topolSlave = getRowsMatrix(surfSlave.connectivity,1:surfSlave.num);


      for vtkSlave = surfSlave.vtkTypes

        elSlave = getElement(obj,vtkSlave,s);

        for vtkMaster = surfSlave.vtkTypes

          elMaster = getElement(obj,vtkSlave,m);

          % loop over pairs of connected master/slave elements
          for iPair = 1:obj.quadrature.numbInterfacePairs

            is = elemPairs(iPair,s);
            im = elemPairs(iPair,m);

            if surfSlave.VTKType(is) ~= vtkSlave; continue; end
            if surfMaster.VTKType(im) ~= vtkMaster; continue; end

            contactState = obj.activeSet.curr(is);

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

            % equilibrium equation and stabilization work with delta traction
            dTrac = trac - tIni;

            nodeMaster = surfMaster.loc2glob(topolMaster(im,1:elMaster.nNode));
            umDof = dofMaster.getLocalDoF(fldM,nodeMaster);

            [Nslave,Nmaster,Nmult] = ...
              getMortarBasisFunctions(obj.quadrature,im,is,elMaster,elSlave,xiMaster,xiSlave);

            % reshape basis function matrices to match number of components
            [Ns,Nm,Nmult] = reshapeBasisFunctions(3,Nslave,Nmaster,Nmult);

            % rotation matrix
            R = getRotationMatrix(obj,MortarSide.slave,is);

            % normal gap (minus sign to be checked)
            g_n = state.gap(3*is-2);

            % tangential slip
            dgt = slip([3*is-1; 3*is]);

            dgtStab = tangSlip([2*is-1; 2*is]);
            slipNorm = norm(dgtStab);

            % operator mapping global vectors to local tangential coordinates
            T = (R(:,2:3))';

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

            % assemble jacobian and rhs of traction balance equations

            % STICK MODE
            if contactState == ContactMode.stick

              asbMt.localAssembly(tDof,umDof,Aum');
              asbDt.localAssembly(tDof,usDof,-Aus');

              % this term can be assembled easily without numerical integration
              rhsT(tDof(1)) = rhsT(tDof(1)) + area*g_n;
              rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + area*dgt;

            end

            % SLIP MODE
            if contactState == ContactMode.slip || contactState == ContactMode.newSlip

              % Independent-variable regularized Coulomb equation:
              %
              %   r_T = t_T - tau_lim*d(t_T,g_T) = 0,
              %
              % with
              %
              %   d = normalize(alpha*n_t + (1-alpha)*n_g),
              %   n_t = normalize(t_T),
              %   n_g = normalize(g_T),
              %   alpha = alpha_min + (1-alpha_min)*(1+n_t'*n_g)/2.
              %
              % No penalty trial traction is used here: t_T and g_T are
              % independent unknowns.

              [tT_lim,dtdgt,dtdtT,dtdtn] = ...
                computeRegularizedFrictionUpdate(obj,trac,dgtStab);

              % Normal contact equation: g_N = 0 in slip.
              asbMt.localAssembly(tDof(1),umDof,Aum(:,1));
              asbDt.localAssembly(tDof(1),usDof,-Aus(:,1));

              % Tangential residual derivative wrt displacement through g_T:
              % d(r_T)/du = -d(tT_lim)/d(g_T) * d(g_T)/du.
              Atu_m = MortarQuadrature.integrate(f2,dtdgt,pagemtimes(T,Nm),dJw);
              Atu_s = MortarQuadrature.integrate(f2,dtdgt,pagemtimes(T,Ns),dJw);
              asbMt.localAssembly(tDof(2:3),umDof,-Atu_m);
              asbDt.localAssembly(tDof(2:3),usDof,Atu_s);

              % Tangential residual derivative wrt normal traction:
              % d(r_T)/d(t_N) = -d(tT_lim)/d(t_N).
              Atn = area*dtdtn;
              asbQ.localAssembly(tDof(2:3),tDof(1),-Atn);

              % Tangential residual derivative wrt tangential traction:
              % d(r_T)/d(t_T) = I - d(tT_lim)/d(t_T).
              Att = area*(eye(2)-dtdtT);
              asbQ.localAssembly(tDof(2:3),tDof(2:3),Att);

              rhsT(tDof(1)) = rhsT(tDof(1)) + area*g_n;

              % Enforce the regularized Coulomb equation.
              rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + area*(trac(2:3)-tT_lim);

            end

            % OPEN MODE
            if contactState == ContactMode.open

              % A_oo
              Aoo = MortarQuadrature.integrate(f1,Nmult,Nmult,dJw);
              asbQ.localAssembly(tDof,tDof,Aoo);

              % rhs (mu,t)
              rhsT(tDof) = rhsT(tDof) + area*dTrac;
            end

          end % end inner master elems loop

        end
      end

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
      nq = ncomp*N2;

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

   function [tT_lim,dtdgt,dtdtT,dtdtn] = computeRegularizedFrictionUpdate(obj,trac,gT)

  % Regularized friction update for independent traction/gap variables.
  %
  % alpha weights the gap direction:
  %
  %   q = alpha*ng + (1-alpha)*nt
  %
  % with
  %
  %   alpha = alphaMin + (alphaMax-alphaMin)*(1-nt'*ng)/2
  %
  % No smooth normalization is used for nt and ng. If either norm is too
  % small, the potentially problematic directional Jacobian terms are
  % skipped.

  tanPhi = tan(deg2rad(obj.phi));

  tN = trac(1);
  tT = trac(2:3);

  tauLim = obj.cohesion - tanPhi*tN;

  tolT = 1e-12;
  tolG = 1e-12;
  tolQ = 1e-12;

  normT = norm(tT);
  normG = norm(gT);

  tT_lim = zeros(2,1);
  dtdgt  = zeros(2,2);
  dtdtT  = zeros(2,2);
  dtdtn  = zeros(2,1);

  if tauLim <= 0
    return
  end

  haveT = normT > tolT;
  haveG = normG > tolG;

  if haveT
    nt = tT/normT;
  else
    nt = zeros(2,1);
  end

  if haveG
    ng = gT/normG;
  else
    ng = zeros(2,1);
  end

  % ------------------------------------------------------------
  % Direction selection
  % ------------------------------------------------------------

  if haveT && haveG

    c = nt.'*ng;

    beta = 0.5*(obj.alphaMax - obj.alphaMin);

    alpha = obj.alphaMin + beta*(1.0 - c);
    alpha = min(obj.alphaMax,max(obj.alphaMin,alpha));

    q = alpha*ng + (1.0-alpha)*nt;

  elseif haveG

    alpha = 1.0;
    beta = 0.0;
    q = ng;

  elseif haveT

    alpha = 0.0;
    beta = 0.0;
    q = nt;

  else

    return

  end

  normQ = norm(q);

  if normQ <= tolQ
    if haveG
      d = ng;
    elseif haveT
      d = nt;
    else
      return
    end
  else
    d = q/normQ;
  end

  tT_lim = tauLim*d;

  % ------------------------------------------------------------
  % Directional Jacobian
  % ------------------------------------------------------------

  if normQ <= tolQ
    % Direction is selected by fallback. Avoid problematic Jacobian.
    dtdtn = -tanPhi*d;
    return
  end

  Pd = eye(2) - d*d.';
  DdDq = Pd/normQ;

  if haveT && haveG

    Pt = eye(2) - nt*nt.';
    Pg = eye(2) - ng*ng.';

    DntDtT = Pt/normT;
    DngDgT = Pg/normG;

    DalphaDtT = -beta*(DntDtT.'*ng).';
    DalphaDgT = -beta*(DngDgT.'*nt).';

    DqDtT = (1.0-alpha)*DntDtT + (ng-nt)*DalphaDtT;
    DqDgT = alpha*DngDgT + (ng-nt)*DalphaDgT;

    DdDtT = DdDq*DqDtT;
    DdDgT = DdDq*DqDgT;

    dtdtT = tauLim*DdDtT;
    dtdgt = tauLim*DdDgT;

  elseif haveG

    Pg = eye(2) - ng*ng.';

    DngDgT = Pg/normG;
    DdDgT = DdDq*DngDgT;

    dtdgt = tauLim*DdDgT;

  elseif haveT

    Pt = eye(2) - nt*nt.';

    DntDtT = Pt/normT;
    DdDtT = DdDq*DntDtT;

    dtdtT = tauLim*DdDtT;

  end

  dtdtn = -tanPhi*d;

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


