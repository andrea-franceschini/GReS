classdef SolidMechanicsContact < MeshTying

  % solid mechanics solver using piece-wise constant multipliers

  properties
    phi               % friction angle in deg
    cohesion          % cohesion
    contactHelper
    activeSet
    NLIter = 0
    stickNodes     % boundary nodes where contact state should stay stick
  end


  methods
    function obj = SolidMechanicsContact(id,domains,inputStruct)

      obj@MeshTying(id,domains,inputStruct);

      if obj.multiplierLocation ~= entityField.surface
        error("Interface Solver %s is not implemented for multipliers ..." + ...
          "located at %s. The only available entityField is %s",...
          class(obj),obj.multiplierLocation,entityField.surface)
      end

    end

    function registerInterface(obj,varargin)

      input = varargin{1};

      input = readInput(struct('Coulomb',[],'ActiveSet',missing),input);

      params = readInput(struct('cohesion',[],'frictionAngle',[]),input.Coulomb);

      obj.cohesion = params.cohesion;
      obj.phi = params.frictionAngle;

      nDofsInterface = getNumbDoF(obj);

      obj.state.traction = zeros(nDofsInterface,1);
      obj.state.deltaTraction = zeros(nDofsInterface,1);
      obj.state.iniTraction = obj.state.traction;

      % the gap in global coordinates
      obj.state.gap = zeros(nDofsInterface,1);

      obj.state.normalGap = zeros(round(1/3*nDofsInterface),1);
      obj.state.tangentialGap = zeros(round(2/3*nDofsInterface),1);
      obj.state.tangentialSlip = zeros(round(2/3*nDofsInterface),1);

      obj.stateOld = obj.state;

      N = obj.grids(MortarSide.slave).surfaces.num;
      initializeActiveSet(obj,N,input.ActiveSet);
      
    end

    function updateState(obj,du)

      % traction update
      actMult = getMultiplierDoF(obj);

      obj.state.traction(actMult) = obj.state.traction(actMult) + du(1:obj.nMult);
      obj.state.deltaTraction = obj.state.traction - obj.stateOld.traction;
      obj.NLIter = obj.NLIter + 1;

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

      if gresLog().getVerbosity > 2
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

      obj.NLIter = 0;

      oldActiveSet = obj.activeSet.curr;
      surfSlave = obj.grids(MortarSide.slave).surfaces;

      for is = 1:numel(obj.activeSet.curr)

        state = obj.activeSet.curr(is);

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
        t = obj.state.traction(id);
        limitTraction = abs(obj.cohesion - tan(deg2rad(obj.phi))*t(1));

        % report traction during activeSet update
        gresLog().log(4,['\n Element %i: traction: %1.4e %1.4e %1.4e   ' ...
          'Limit tangential traction: %1.4e \n'],is,t(:), limitTraction)

        obj.activeSet.curr(is) = updateContactState(state,t,...
                                                    limitTraction, ...
                                                    obj.state.normalGap(is),...
                                                    obj.activeSet.tol);

      end

      % check if active set changed
      asNew = ContactMode.integer(obj.activeSet.curr);
      asOld = ContactMode.integer(oldActiveSet);

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

      if gresLog().getVerbosity > 2
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

        fprintf('Stick dofs: %i    Slip dofs: %i    Open dofs: %i \n',...
          sum(asNew==1), sum(any([asNew==2,asNew==3],2)), sum(asNew==4));
      end

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
      obj.state.traction = obj.state.traction + tIni;
      obj.state.initTraction = obj.state.traction;

      setStickNodes(obj);

    end


    function addInitialTraction(obj)
      % add a traction on the fault
      obj.state.traction = obj.state.traction + t;
      obj.state.iniTraction = obj.state.iniTraction + t;
    end

    function t = computeInitialTraction(obj)
      % initialize traction for cell stress (average)
      sl = MortarSide.slave;
      poro = obj.domains(sl).getPhysicsSolver("Poromechanics");
      surf = obj.grids(sl).surfaces;
      faces = obj.domains(sl).grid.faces;
      normals = surf.normal;
      faceIds = surf.faceId;
      cellIds = faces.neighbors(faceIds,1);
      sigma = zeros(3);
      idx = [1;6;5;6;2;4;5;4;3];
      t = zeros(getNumbDoF(obj),1);
      for i = 1:numel(cellIds)
        sigma(:) = poro.avStress(i,idx);
        n = normals(i,:);
        tDof = getMultiplierDoF(obj,i);
        t(tDof) = sigma*n';
      end
    end



    function advanceState(obj)

      obj.state.deltaTraction(:) = 0;
      obj.stateOld = obj.state;
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
      obj.activeSet.curr = obj.activeSet.prev;
      obj.state = obj.stateOld;
      obj.state.deltaTraction(:) = 0;
      obj.NLIter = 0;
      if obj.activeSet.resetActiveSet
        resetConfiguration(obj);
      end
    end



    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      outTraction = fac*obj.state.traction + ...
        (1-fac)*obj.stateOld.traction;

      outNormalGap = fac*obj.state.normalGap + ...
        (1-fac)*obj.stateOld.normalGap;
      outTangentialSlip = fac*obj.state.tangentialSlip + ...
        (1-fac)*obj.stateOld.tangentialSlip;
      outTangentialSlip = (reshape(outTangentialSlip,2,[]))';

      outTangentialGap = fac*obj.state.tangentialGap + ...
        (1-fac)*obj.stateOld.tangentialGap;
      outTangentialGap = (reshape(outTangentialGap,2,[]))';

      outTangentialGapNorm = sqrt(outTangentialGap(:,1).^2 + ...
        outTangentialGap(:,2).^2);

      tT = [outTraction(2:3:end),outTraction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      fractureState = ContactMode.integer(obj.activeSet.curr);

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
        };

      surfaceStr = cell2struct(entries, {'name','data'}, 2);
    end

    function writeSolution(obj,fac,tID)

      outTraction = fac*obj.state.traction + ...
        (1-fac)*obj.stateOld.traction;
      outNormalGap = fac*obj.state.normalGap + ...
        (1-fac)*obj.stateOld.normalGap;
      outSlip = fac*obj.state.tangentialSlip + ...
        (1-fac)*obj.stateOld.tangentialSlip;
      outSliding = fac*obj.state.tangentialGap + ...
        (1-fac)*obj.stateOld.tangentialGap;

      tT = [outTraction(2:3:end),outTraction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      obj.outstate.results(tID).tractionVec = outTraction;
      obj.outstate.results(tID).normalGap = outNormalGap;
      obj.outstate.results(tID).slipIncrement = outSlip;
      obj.outstate.results(tID).tangentialGap = outSliding;
      obj.outstate.results(tID).tangentialTractionNorm = norm_tT;

    end

  end

  methods (Access = protected)


    function computeGap(obj)
      % compute normal gap and tangential slip (local coordinates)

      um = obj.domains(MortarSide.master).state.data.displacements;
      us = obj.domains(MortarSide.slave).state.data.displacements;

      % recover variationally consistent stabilized gaps
      areaSlave = repelem(obj.getSlaveArea(),3);

      areaGap = (obj.D*us + obj.M*um);

      obj.state.gap = areaGap./areaSlave;

      [~,rhsStab] = getStabilizationMatrixAndRhs(obj);

      stabGap = (areaGap + rhsStab)./areaSlave;
      stabSlip = (obj.state.gap-obj.stateOld.gap) + rhsStab./areaSlave;

      stabSlip(1:3:end) = [];

      obj.state.tangentialSlip = stabSlip;
      %
      obj.state.normalGap = stabGap(1:3:end);
      obj.state.tangentialGap = obj.stateOld.tangentialGap + stabSlip;

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
      slip = obj.state.gap - obj.stateOld.gap;

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
            trac = obj.state.traction(tDof);

            % equilibrium equation and stabilization work with delta traction
            dTrac = trac - obj.state.iniTraction(tDof);

            nodeMaster = surfMaster.loc2glob(topolMaster(im,1:elMaster.nNode));
            umDof = dofMaster.getLocalDoF(fldM,nodeMaster);

            [Nslave,Nmaster,Nmult] = ...
              getMortarBasisFunctions(obj.quadrature,im,is,elMaster,elSlave,xiMaster,xiSlave);

            % reshape basis function matrices to match number of components
            [Ns,Nm,Nmult] = reshapeBasisFunctions(3,Nslave,Nmaster,Nmult);

            % rotation matrix
            R = getRotationMatrix(obj,MortarSide.slave,is);

            % normal gap (minus sign to be checked)
            g_n = obj.state.gap(3*is-2);

            % tangential slip
            dgt = slip([3*is-1 3*is]);
            slipNorm = norm(dgt);

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

              slidingTol = obj.activeSet.tol.sliding;

              tauLim = obj.cohesion - trac(1)*tan(deg2rad(obj.phi));

              asbMt.localAssembly(tDof(1),umDof,Aum(:,1));
              asbDt.localAssembly(tDof(1),usDof,-Aus(:,1));

              % A_tu (non linear term)
              if slipNorm > slidingTol && obj.NLIter > 0

                % compute only on slip terms with sliding large enough
                dtdgt = computeDerTracGap(obj,trac(1),dgt);
                Atu_m = MortarQuadrature.integrate(f2, dtdgt,pagemtimes(T,Nm),dJw);
                Atu_s = MortarQuadrature.integrate(f2, dtdgt,pagemtimes(T,Ns),dJw);
                asbMt.localAssembly(tDof(2:3),umDof,-Atu_m);
                asbDt.localAssembly(tDof(2:3),usDof,Atu_s);

                % A_tn (non linear term)
                dtdtn = computeDerTracTn(obj,dgt);
                Atn = area*dtdtn;
                asbQ.localAssembly(tDof(2:3),tDof(1),-Atn);

                tT_lim = tauLim*(dgt/norm(dgt));

              else

                % if slip is small, use current traction
                vaux = trac(2:3);
                dtdtn = - tan(deg2rad(obj.phi))*vaux/norm(vaux);
                Atn = area*dtdtn;
                asbQ.localAssembly(tDof(2:3),tDof(1),-Atn);

                tT_lim = tauLim*vaux/norm(vaux);

              end

              % A_tt
              Att = area*eye(2);
              asbQ.localAssembly(tDof(2:3),tDof(2:3),Att);

              rhsT(tDof(1)) = rhsT(tDof(1)) + area*g_n;

              % rhs (mu_t,tT) - local frame
              rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + area * (trac(2:3)-tT_lim);


              if gresLog().getVerbosity > 5
                fprintf('\nelement %i- rhsT: %5.3e %5.3e \n',is,Att*trac(2:3))
                fprintf('element %i- rhsTlim: %5.3e %5.3e \n',is,MortarQuadrature.integrate(f1,Nmult_t,tT_lim,dJw))
                fprintf('------------------------------------ \n')
              end

            end

            % OPEN MODE
            if contactState == ContactMode.open

              % A_oo
              Aoo = MortarQuadrature.integrate(f1,Nmult,Nmult,dJw);
              asbQ.localAssembly(tDof,tDof,Aoo);

              % rhs (mu,t)
              rhsT(tDof) = rhsT(tDof) + area*trac;
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
      rhsH = -H*obj.state.deltaTraction;

      rhsH(1:3:end) = -H(1:3:end,:) * obj.state.traction;

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

    function dtdgt = computeDerTracGap(obj,sigma_n,slip)
      % gt = obj.g_T(get_dof(nodeId));

      % slip: 2x1 local tangential slip
      tauLim = obj.cohesion - tan(deg2rad(obj.phi))*sigma_n;
      dtdgt = tauLim*((eye(2)*norm(slip)^2 - slip*slip')/(norm(slip))^3);

    end

    function dtdtn = computeDerTracTn(obj,slip)

      tanPhi = tan(deg2rad(obj.phi));

      if norm(slip) > obj.activeSet.tol.sliding
        %use available gap to properly compute traction
        dtdtn = -tanPhi*(slip/norm(slip));
      else
        t = obj.state.traction(getMultiplierDoF(obj,is));
        t = t(2:3);
        % convert to global reference system
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


