classdef SolidMechanicsContact < MeshTying

  % solid mechanics solver using piece-wise constant multipliers

  properties
    phi               % friction angle in deg
    cohesion          % cohesion
    contactHelper
    activeSet
    oldStab = false
    NLIter = 0
  end


  methods
    function obj = SolidMechanicsContact(id,domains,inputStruct)

      obj@MeshTying(id,domains,inputStruct);

    end

    function registerInterface(obj,varargin)

      input = varargin{1};

      obj.stabilizationScale = getXMLData(input.Quadrature,1.0,"stabilizationScale");
      obj.cohesion = getXMLData(input.Coulomb,[],"cohesion");
      obj.phi = getXMLData(input.Coulomb,[],"frictionAngle");

 
      nDofsInterface = getNumbDoF(obj);

      obj.state.traction = zeros(nDofsInterface,1);
      obj.state.deltaTraction = zeros(nDofsInterface,1);
      obj.state.iniTraction = obj.state.traction;

      % the gap in global coordinates
      obj.state.gap = zeros(nDofsInterface,1);

      obj.state.normalGap = zeros(round(1/3*nDofsInterface),1);
      obj.state.tangentialGap = zeros(round(1/3*nDofsInterface),1);
      obj.state.tangentialSlip = zeros(round(2/3*nDofsInterface),1);

      obj.stateOld = obj.state;

      N = getMesh(obj,MortarSide.slave).nSurfaces;
      SolidMechanicsContact.initializeActiveSet(obj,input,N);

    end

    function updateState(obj,du)

      % traction update
      actMult = getMultiplierDoF(obj);
      %obj.state.deltaTraction(actMult) = du(1:obj.nMult);
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
      % if ~obj.oldStab
      %   computeStabilizationMatrix(obj);
      % end

      % get stabilization matrix depending on the current active set
      [H,rhsStab] = getStabilizationMatrixAndRhs(obj);
      %[H,rhsStab] = getStabilizationMatrixAndRhsOld(obj);
      

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



    function hasChanged = updateActiveSet(obj)

      obj.NLIter = 0;

      oldActiveSet = obj.activeSet.curr;
      mshSlave = getMesh(obj,MortarSide.slave);

      for is = 1:numel(obj.activeSet.curr)

        state = obj.activeSet.curr(is);

        % force boundary element to stick state
        nodes = obj.interfMesh.local2glob{2}(mshSlave.surfaces(is,:));

        if obj.activeSet.forceStickBoundary
          if any(ismember(nodes,obj.dirNodes))
            % element has a dirichlet node - keep it stick
            continue
          end
        end

        id = DoFManager.dofExpand(is,3);
        t = obj.state.traction(id);
        limitTraction = abs(obj.cohesion - tan(deg2rad(obj.phi))*t(1));

        % report traction during activeSet update
        gresLog().log(3,['\n Element %i: traction: %1.4e %1.4e %1.4e   ' ...
          'Limit tangential traction: %1.4e \n'],is,t(:), limitTraction)

        obj.activeSet.curr(is) = obj.updateContactState(state,t,...
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
      

      hasChanged = any(diffState);

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

      if hasChanged


        % EXCEPTION 1): check if area of fracture changing state is relatively small
        msh = getMesh(obj,MortarSide.slave);
        areaChanged = sum(msh.surfaceArea(hasChangedElem));
        totArea = sum(msh.surfaceArea);
        if areaChanged/totArea < obj.activeSet.tol.areaTol
          %obj.activeSet.curr = oldActiveSet;
          % change the active set, but flag it as nothing changed
          hasChanged = false;
          gresLog().log(1,['Active set update suppressed due to small fracture change:' ...
            ' areaChange/areaTot = %3.2e \n'],areaChanged/totArea);
        end

        % EXCEPTION 2): check if changing elements have been looping from
        % stick to slip/open too much times

        if all(obj.activeSet.stateChange(hasChangedElem) > obj.activeSet.tol.maxStateChange)
          hasChanged = false;
          gresLog().log(1,['Active set update suppressed due to' ...
            ' unstable behavior detected'])
        end
      end

      % if hasChanged && ~obj.oldStab
      computeStabilizationMatrix(obj);
      % end
    end



    function advanceState(obj)

      obj.state.deltaTraction(:) = 0;
      obj.stateOld = obj.state;
      obj.activeSet.prev = obj.activeSet.curr;
      obj.NLIter = 0;

      % reset the counter for changed states
      obj.activeSet.stateChange(:) = 0;


      %resetConfiguration(obj);

      % if obj.resetActiveSet
      %   % reset everything to stick for next time step
      %   obj.activeSet.curr(:) = ContactMode.stick;
      % end
    end

    function isReset = resetConfiguration(obj)
      if ~obj.activeSet.resetActiveSet
        isReset = false;
      else
        % obj.state = obj.stateOld;
        toReset = obj.activeSet.curr(:) ~= ContactMode.open;
        obj.activeSet.curr(toReset) = ContactMode.stick;

        % recompute stabilization matrix
        computeStabilizationMatrix(obj);
        isReset = true;
      end
    end



    function goBackState(obj,dt)
      % reset state to beginning of time step
      %updateActiveSet(obj);
      obj.activeSet.curr = obj.activeSet.prev;
      obj.state = obj.stateOld;
      obj.state.deltaTraction(:) = 0;
      obj.NLIter = 0;
      %resetConfiguration(obj);
    end



    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      outTraction = fac*obj.state.traction + ...
        (1-fac)*obj.stateOld.traction;


      outNormalGap = fac*obj.state.normalGap + ...
        (1-fac)*obj.stateOld.normalGap;
      outTangentialSlip = fac*obj.state.tangentialSlip + ...
        (1-fac)*obj.stateOld.tangentialSlip;
      outTangentialSlip = (reshape(outTangentialSlip,2,[]))';

      outSliding = fac*obj.state.tangentialGap + ...
        (1-fac)*obj.stateOld.tangentialGap;

      outSlipNorm = sqrt(outTangentialSlip(:,1).^2 + outTangentialSlip(:,2).^2);

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
        'tangential_slip_norm',    outSlipNorm
        'tangential_gap_norm',     outSliding
        'fracture_state',          fractureState
        };

      surfaceStr = cell2struct(entries, {'name','data'}, 2);
    end

    function writeMatFile(obj,fac,tID)

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

      %nS = obj.interfMesh.msh(2).nSurfaces;

      um = obj.domains(1).state.data.displacements;
      us = obj.domains(2).state.data.displacements;
      % 


      % recover variationally consistent stabilized gaps
      areaSlave = repelem(obj.getSlaveArea(),3);

      areaGap = (obj.D*us + obj.M*um);

      obj.state.gap = areaGap./areaSlave;

      %[~,rhsStab] = getStabilizationMatrixAndRhs(obj);
      %stabGap = -rhsH./areaSlave;
      %gap = H*(obj.state.traction - obj.state.iniTraction)./areaSlave;

      %gap = obj.state.gap - stabGap;
      slip = obj.state.gap - obj.stateOld.gap;
      slip(1:3:end) = [];

      obj.state.tangentialSlip = slip;
      % 
      obj.state.normalGap = obj.state.gap(1:3:end);
      obj.state.tangentialGap = sqrt(obj.state.gap(2:3:end).^2 + ...
        obj.state.gap(3:3:end).^2);

      % rotate gap from global to local coordinates
      %       for i = 1:nS
      %         R  = getRotationMatrix(obj,i);
      %         obj.dispJump.curr(dofId(i,3)) = R'*gap(dofId(i,3));
      %         locSlip = R(:,2:3)'*slipIncrement(dofId(i,3));
      %         obj.slip.curr(i) = norm(locSlip);
      %       end


      %       if obj.domains(2).simparams.verbosity > 2
      %         m = max(abs(stabGap./(obj.Jst*us)));
      %         fprintf('Maximum gap deviation due to stabilization: %4.1f %% \n',100*m)
      %       end

      %       gapOld = (obj.D*us_old - obj.M*um_old - stabGap)./sum(obj.D,2);

    end


    function computeContactMatricesAndRhs(obj)

      % Compute contact matrices and rhs

      % the syntax of the comments is inspired by the appendix of
      % Franceschini, Andrea, et al. "Algebraically stabilized Lagrange
      % multiplier method for frictional contact mechanics with
      % hydraulically active fractures." Computer Methods in Applied
      % Mechanics and Engineering 368 (2020): 113161.

      dofSlave = getDoFManager(obj,MortarSide.slave);
      dofMaster = getDoFManager(obj,MortarSide.master);

      % extract current and previous displacements
      % us = obj.domains(2).state.data.displacements;
      % us_old = obj.domains(2).stateOld.data.displacements;
      % um = obj.domains(1).state.data.displacements;
      % um_old = obj.domains(1).stateOld.data.displacements;

      % compute displacement increment
      %       um_slip = um - um_old;
      %       us_slip = us - us_old;

      % define matrix assemblers
      [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj);

      % define rhs vectors;
      rhsUm = zeros(getNumbDoF(dofMaster),1);
      rhsUs = zeros(getNumbDoF(dofSlave),1);
      rhsT = zeros(getNumbDoF(obj),1);

      fldM = dofMaster.getVariableId(obj.coupledVariables);
      fldS = dofSlave.getVariableId(obj.coupledVariables);

      % anonymous functions for local fem computations
      f1 = @(a,b) pagemtimes(a,'ctranspose',b,'none');
      f2 = @(a,b,c) pagemtimes(a,'transpose',pagemtimes(b,c),'none');
      f3 = @(a,b) pagemtimes(a,b);

      % compute slip
      slip = obj.state.gap - obj.stateOld.gap;


      for iPair = 1:obj.quadrature.numbInterfacePairs

        is = obj.quadrature.interfacePairs(iPair,1);
        im = obj.quadrature.interfacePairs(iPair,2);

        contactState = obj.activeSet.curr(is);

        % retrieve mortar integration data
        xiMaster = obj.quadrature.getMasterGPCoords(iPair);
        xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
        dJw = obj.quadrature.getIntegrationWeights(iPair);

        area = sum(dJw);

        % define slave related quantities
        nodeSlave = obj.interfMesh.local2glob{2}(obj.interfMesh.msh(2).surfaces(is,:));
        usDof = dofSlave.getLocalDoF(fldS,nodeSlave);
        tDof = getMultiplierDoF(obj,is);
        trac = obj.state.traction(tDof);
        % equilibrium equation and stabilization work with delta traction
        dTrac = trac - obj.state.iniTraction(tDof);

        % displacement increment to check if this is the first newton
        % iteration
        % usDiff = norm(us(usDof) - us_old(usDof));

        nodeMaster = obj.interfMesh.local2glob{1}(obj.interfMesh.msh(1).surfaces(im,:));
        umDof = dofMaster.getLocalDoF(fldM,nodeMaster);

        [Nslave,Nmaster,Nmult] = ...
          getMortarBasisFunctions(obj.quadrature,im,is,xiMaster,xiSlave);

        % reshape basis function matrices to match number of components
        [Ns,Nm,Nmult] = reshapeBasisFunctions(3,Nslave,Nmaster,Nmult);

        % Nmult(1,1,:) = -Nmult(1,1,:); 

        % rotation of multipliers
        R = getRotationMatrix(obj.interfMesh,is);
        %NmultR = pagemtimes(Nmult,R);

        % Reduce the dimension of multiplier basis functions exploiting
        % the local definition of degrees of freedom

        % the normal component of the multipliers basis
        %Nmult_n = -Nmult(1,1,:);

        % tangential component of multiplier basis
        %Nmult_t = NmultR(:,2:3,:);

        % get normal of the slave element
        %normalEl = obj.interfMesh.getNormal(is);
        
        %Tel = eye(3) - normalEl*normalEl';

        % get normal at the gauss points (for warped facets...?)
        % normalNodes = obj.interfMesh.getNodalNormal(is);
        % normal = pagemtimes(Ns,normalNodes);

        % operator selecting only tangential components of the
        % displacements
        %T = eye(3) - pagemtimes(normal,'none',normal,'transpose');
        % normal and tangential component of displacement basis functions
        % Nm_n = pagemtimes(normal,'transpose',Nm,'none');
        % Nm_t = pagemtimes(T,Nm);
        % 
        % Ns_n = pagemtimes(normal,'transpose',Ns,'none');
        % Ns_t = pagemtimes(T,Ns);


        % normal gap (minus sign to be checked)
        g_n = obj.state.normalGap(is);

        % tangential slip (quasi-static coulomb law) in local coords
        %tangDofs1 = [2*is-1 2*is];

        % operator mapping global vectors to local tangential coordinates
        T = (R(:,2:3))';
        dgt = slip([3*is-1 3*is]);
        slipNorm = norm(dgt);



        % slip is stabilized while does not incorporate the stabilizing
        % contribution, hence can be non zero even for new slip elements

        % A_us
        Aum =  MortarQuadrature.integrate(f1,Nm,Nmult,dJw);
        Aus =  MortarQuadrature.integrate(f1,Ns,Nmult,dJw);
        % apply rotation matrix due to mixed dof assembly
        Aum = Aum*R;
        Aus = Aus*R;
        asbMu.localAssembly(umDof,tDof,Aum);
        asbDu.localAssembly(usDof,tDof,-Aus);

        % Aut_m =  MortarQuadrature.integrate(f1,Nm_t,Nmult_t,dJw);
        % Aut_s =  MortarQuadrature.integrate(f1,Ns_t,Nmult_t,dJw);
        % asbMu.localAssembly(umDof,tDof(2:3),-Aut_m);
        % asbDu.localAssembly(usDof,tDof(2:3),Aut_s);

        % rhs (jump(eta),t)
        rhsUm(umDof) = rhsUm(umDof) + Aum*dTrac;
        rhsUs(usDof) = rhsUs(usDof) - Aus*dTrac;

        % assemble jacobian and rhs of traction balance equations

        % STICK MODE
        if contactState == ContactMode.stick

          asbMt.localAssembly(tDof,umDof,Aum');
          asbDt.localAssembly(tDof,usDof,-Aus');

          % this term can be assembled easily without numerical integration
          % (mu_n,g_n)
          rhsT(tDof(1)) = rhsT(tDof(1)) + area*g_n;
          rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + area*dgt;

        end

        % SLIP MODE
        if contactState == ContactMode.slip || contactState == ContactMode.newSlip

          slidingTol = obj.activeSet.tol.sliding;

          % A_nu
          % Anu_m = MortarQuadrature.integrate(f1, Nmult_n, Nm_n, dJw);
          % Anu_s = MortarQuadrature.integrate(f1, Nmult_n, Ns_n, dJw);
          asbMt.localAssembly(tDof(1),umDof,Aum(:,1));
          asbDt.localAssembly(tDof(1),usDof,-Aus(:,1));

          % A_tu (non linear term)
          if slipNorm > slidingTol && obj.NLIter > 0 
            % compute only on slip terms with sliding large enough
            dtdgt = computeDerTracGap(obj,trac(1),dgt);
            Atu_m = MortarQuadrature.integrate(f3, dtdgt,pagemtimes(T,Nm),dJw);
            Atu_s = MortarQuadrature.integrate(f3, dtdgt,pagemtimes(T,Ns),dJw);
            asbMt.localAssembly(tDof(2:3),umDof,-Atu_m);
            asbDt.localAssembly(tDof(2:3),usDof,Atu_s);

            % A_tn (non linear term)
            dtdtn = computeDerTracTn(obj,dgt);
            Atn = area*dtdtn; 
            asbQ.localAssembly(tDof(2:3),tDof(1),-Atn);
          end

          % A_tt
          Att = area*eye(2);
          asbQ.localAssembly(tDof(2:3),tDof(2:3),Att);


          % rhs (mu_n,g_n)
          % rhsT(tDof(1)) = rhsT(tDof(1)) + ...
          %   MortarQuadrature.integrate(f1,Nmult_n,g_n,dJw);
          rhsT(tDof(1)) = rhsT(tDof(1)) + area*g_n;

          % rhs -(mu_t,t*_T) (non linear term)
          tT_lim = computeLimitTraction(obj,dgt,trac,slipNorm);


          % rhs (mu_t,tT) tT and tT_lim in local coordinates (to be consistent with dof definition)
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

      % assemble matrices into jacobian blocks
      obj.addJum(MortarSide.master, asbMu.sparseAssembly());
      obj.addJum(MortarSide.slave, asbDu.sparseAssembly());
      obj.addJmu(MortarSide.master, asbMt.sparseAssembly());
      obj.addJmu(MortarSide.slave, asbDt.sparseAssembly());

      obj.Jconstraint = asbQ.sparseAssembly();

      obj.addRhs(MortarSide.master,rhsUm);
      obj.addRhs(MortarSide.slave,rhsUs);
      obj.rhsConstraint = rhsT;

      % assemble rhs contribution

      % obj.Jmu = asbMu.sparseAssembly();
      % obj.Jsu = asbDu.sparseAssembly();
      % obj.Jmt = asbMt.sparseAssembly();
      % obj.Jst = asbDt.sparseAssembly();
      % obj.Jtt = asbQ.sparseAssembly();

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
      if obj.oldStab
        elOpen = find(obj.activeSet.curr == ContactMode.open);
        elSlip = [find(obj.activeSet.curr == ContactMode.slip);...
          find(obj.activeSet.curr == ContactMode.newSlip)];

        dofOpen = DoFManager.dofExpand(elOpen,3);
        dofSlip = [3*elSlip-1; 3*elSlip];

        % remove rows and columns of dofs not requiring stabilization
        H([dofOpen;dofSlip],:) = 0;
        H(:,[dofOpen;dofSlip]) = 0;

        % scale the off-diagonal terms of the stabilization matrix to satisfy local
        % conservation
        % dH = diag(H);
        % H = H - diag(dH);
        % f = abs(dH./sum(H,2));
        % f(isnan(f) | isinf(f)) = 0; 
        % H = f.*H + diag(dH);

        % rowSum = sum(H,2);
        % H = H - diag(rowSum);
        % 
        % if ~obj.oldStab
        %   assert(norm(sum(H,2))<1e-8, 'Stabilization matrix is not locally conservative')
        % end

      end
      rhsH = -H*obj.state.deltaTraction;
      rhsH(1:3:end) = -H(1:3:end,:) * obj.state.traction;

      % check that jump stabilizing terms are properly removed
      if ~obj.oldStab
        assert(norm(sum(H,2))<1e-8, 'Stabilization matrix is not locally conservative')
      end

    end

    function computeStabilizationMatrix(obj)


      nComp = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);

      % initialize matrix estimating number of entries
      % number of internal slave elements
      nes = sum(all(obj.interfMesh.e2f{2},2));
      nEntries = 2*36*nes; % each cell should contribute at least two times
      nmult = getNumbDoF(obj);
      localKernel = @(S,e1,e2) assembleLocalStabilization(obj,S,e1,e2);
      asbH = assembler(nEntries,nmult,nmult,localKernel);

      % get list of internal master edges
      inEdgeMaster = find(all(obj.interfMesh.e2f{1},2));

      for ieM = inEdgeMaster'
        % loop over internal master edges

        % get 2 master faces sharing internal edge ie
        fM = obj.interfMesh.e2f{1}(ieM,:);
        assert(numel(fM)==2,['Unexpected number of connected faces for' ...
          'master edge %i. Expected 2.'], ieM);

        % get slave faces sharing support with 2 master faces
        fS = unique([find(obj.interfMesh.elemConnectivity(fM(1),:)),...
          find(obj.interfMesh.elemConnectivity(fM(2),:))]);

        if numel(fS) < 2
          continue
        end

        % average master elements area
        Am = mean(obj.interfMesh.msh(1).surfaceArea(fM));

        % get internal edges of slave faces
        eS = unique(obj.interfMesh.f2e{2}(fS,:));
        id = all(ismember(obj.interfMesh.e2f{2}(eS,:),fS),2);
        ieS = eS(id);

        % get master/slave nodes in the local macroelement
        nM = obj.interfMesh.e2n{1}(ieM,:);
        nS = unique(obj.interfMesh.e2n{2}(eS,:));

        % compute local schur complement approximation
        S = computeSchurLocal(obj,nM,nS,fS);

        % assemble stabilization matrix component
        for iesLoc = ieS'

          % loop over internal slave edges in the macroelement

          % get pair of slave faces sharing the edge
          f = obj.interfMesh.e2f{2}(iesLoc,:);

          % mean area of the slave faces
          As = mean(obj.interfMesh.msh(2).surfaceArea(f));
          fLoc1 = dofId(find(fS==f(1)),nComp);
          fLoc2 = dofId(find(fS==f(2)),nComp);

          % local schur complement for macroelement pair of slave faces
          Sloc = 0.5*(Am/As)*(S(fLoc1,fLoc1)+S(fLoc2,fLoc2));
          asbH.localAssembly(Sloc,f(1),f(2));
        end
      end

      obj.stabilizationMat = asbH.sparseAssembly();

      %assert(norm(sum(obj.stabilizationMat,2))<1e-8, 'Stabilization matrix is not locally conservative')
    end

    function [dofRow,dofCol,mat] = assembleLocalStabilization(obj,S,e1,e2)
      % assemble stabilization matrix S (in global coordinates) for
      % elements e1 and e2.

      % for contact problem, we need to carefully handle element in mixed
      % contact state

      contactState = obj.activeSet.curr([e1,e2]);

      nc = getDoFManager(obj,MortarSide.slave).getNumberOfComponents(obj.coupledVariables);
      dof1 = DoFManager.dofExpand(e1,nc);
      dof2 = DoFManager.dofExpand(e2,nc);

      if nc > 1
        % vector field, rotation matrix needed

        % get average rotation matrix
        n1 = getNormal(obj.interfMesh,e1);
        n2 = getNormal(obj.interfMesh,e2);
        if abs(n1'*n2 -1) < 1e4*eps
          avgR = obj.interfMesh.computeRot(n1);
        else
          A1 = obj.interfMesh.msh(2).surfaceArea(e1);
          A2 = obj.interfMesh.msh(2).surfaceArea(e2);
          nAvg = n1*A1 + n2*A2;
          nAvg = nAvg/norm(nAvg);
          avgR = obj.interfMesh.computeRot(nAvg);
        end

        % apply rotation matrix to S
        %S = avgR'*S*avgR;

      end

      % prepare matrix for full stick edge
      mat = obj.stabilizationScale*[S,-S;-S,S];
      dofRow = [dof1;dof2];
      dofCol = [dof1;dof2];

      if any(contactState == ContactMode.open)
        % no stabilization required for open dofs
        mat(:) = 0;
      elseif any(contactState ~= ContactMode.stick)
        % mixed state elements
        mat([2,3,5,6],:) = 0;
        mat(:,[2,3,5,6]) = 0;
        % if contactState(1) ~= ContactMode.stick
        %   mat([2 3],:) = 0;
        %   %mat(:,[2 3]) = 0;
        % end
        % if contactState(2) ~= ContactMode.stick
        %   mat([5 6],:) = 0;
        %   %mat(:,[5 6]) = 0;
        % end
      end

    end


    function [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj)
      % helper to define contact matrix assemblers

      dofSlave = getDoFManager(obj,MortarSide.slave);
      dofMaster = getDoFManager(obj,MortarSide.master);

      ncomp = 3;

      % get number of index entries for sparse matrices

      % number of master nodes attacched to slave elements
      nNmaster = obj.interfMesh.msh(1).surfaceNumVerts'*obj.interfMesh.elemConnectivity;

      N1 = sum(nNmaster);
      N2 = sum(obj.interfMesh.elemConnectivity,1)*obj.interfMesh.msh(2).surfaceNumVerts;

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

      % return the limit traction vector in the global frame
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
  end


  methods (Static)

    function var = getCoupledVariables()
      var = Poromechanics.getField();
    end


    function initializeActiveSet(contactSolver,input,N)
      % initialize an Active Set of size N and write it on input obj
      % read xml field ActiveSet
      
      assert(isprop(contactSolver,"activeSet"),"A property named 'activeSet' " + ...
        "is required to initialize the activeSet");

      contactSolver.activeSet.curr = repmat(ContactMode.stick,N,1);
      contactSolver.activeSet.prev = contactSolver.activeSet.curr;
      % count how many times an element has changed state during iteration
      contactSolver.activeSet.stateChange = zeros(N,1);

      % extract xml field
      if isfield(input,"ActiveSet")
        input = input.ActiveSet;
      end

      % optional flags
      contactSolver.activeSet.resetActiveSet = ...
        logical(getXMLData(input,1,"resetActiveSet"));
      contactSolver.activeSet.forceStickBoundary = ...
        logical(getXMLData(input,0,"forceStickBoundary"));

      % tolerances
      if isfield(input,"Tolerances")
        input = input.Tolerances;
      end
      contactSolver.activeSet.tol.sliding = getXMLData(input,1e-4,"sliding");
      contactSolver.activeSet.tol.normalGap = getXMLData(input,1e-6,"normalGap");
      contactSolver.activeSet.tol.normalTrac = getXMLData(input,1e-3,"normalTraction");
      contactSolver.activeSet.tol.slidingCheck = getXMLData(input,1e-4,"tangentialViolation");
      contactSolver.activeSet.tol.minLimitTraction = getXMLData(input,1e-4,"minLimitTraction");
      contactSolver.activeSet.tol.areaTol = getXMLData(input,1e-2,"areaChange");
      contactSolver.activeSet.tol.maxStateChange = getXMLData(input,100,"maxActiveSetChange");

    end

    function outState = updateContactState(inState,traction,tLimit,normalGap,tols)
      % tols: structure with tolerances named
      % according to initializeActiveSet method
      % the method is static to be reused by other contact solvERS
      outState = inState;


      % contact state update
      if inState == ContactMode.open
        % from open to stick
        if normalGap < -  tols.normalGap
          outState = ContactMode.stick;
        end

      elseif traction(1) > tols.normalTrac
        outState = ContactMode.open;

      else % not open

        % tangential traction norm
        tau = norm(traction(2:3));
        % limiting traction

        % set to 0 a tLimit that is too small
        if tLimit < tols.minLimitTraction
          tLimit = 0;
        end

        % relax stick/sliding transition if the state is about to change
        if inState == ContactMode.stick && tau >= tLimit

          % reduce the tau if goes above limit
          tau = tau*(1-tols.slidingCheck);

        elseif inState ~= ContactMode.stick  && tau <=tLimit

          % increase tau if falls below limit
          tau = tau*(1+tols.slidingCheck);
        end

        % change the state after relaxation
        if tau > tLimit
          if inState == ContactMode.stick
            outState = ContactMode.newSlip;
          else
            outState = ContactMode.slip;
          end
        else
          outState = ContactMode.stick;
        end
      end
    end

  end

end


