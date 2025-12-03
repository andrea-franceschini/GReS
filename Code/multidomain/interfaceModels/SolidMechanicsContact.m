classdef SolidMechanicsContact < MeshTying

  % solid mechanics solver using piece-wise constant multipliers

  properties
    phi               % friction angle in deg
    cohesion          % cohesion
    contactHelper
    activeSet
  end

  properties (Access=protected)
    resetActiveSet = false
    forceStickBoundary = true
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
      obj.state.iniTraction = obj.state.traction;
      obj.state.slip = zeros(round(1/3*nDofsInterface),1);
      obj.state.normalGap = zeros(round(1/3*nDofsInterface),1);
      obj.state.tangentialGap = zeros(round(1/3*nDofsInterface),1);

      obj.stateOld = obj.state;

      N = getMesh(obj,MortarSide.slave).nSurfaces;
      SolidMechanicsContact.initializeActiveSet(obj,input,N);

    end

    function updateState(obj,du)

      % traction update
      actMult = getMultiplierDoF(obj);
      obj.state.traction(actMult) = obj.state.traction(actMult) + du(1:obj.nMult);

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

      computeStabilizationMatrix(obj)

      % get stabilization matrix depending on the current active set
      [H,rhsStab] = getStabilizationMatrixAndRhs(obj);

      obj.Jconstraint = obj.Jconstraint - H;
      obj.rhsConstraint = obj.rhsConstraint + rhsStab;

      if gresLog().getVerbosity > 1
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

      oldActiveSet = obj.activeSet.curr;

      for is = 1:numel(obj.activeSet.curr)

        state = obj.activeSet.curr(is);

        % force boundary element to stick state
        nodes = obj.interfMesh.msh(2).surfaces(is,:);

        if obj.forceStickBoundary
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

      diffState = asNew - asOld;

      idNewSlipToSlip = all([asOld==2 diffState==1],2);
      diffState(idNewSlipToSlip) = 0;
      hasChangedElem = diffState~=0;
      obj.activeSet.stateChange(hasChangedElem) = ...
        obj.activeSet.stateChange(hasChangedElem) + 1;

      hasChanged = any(diffState);

      if gresLog().getVerbosity > 1
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
        % check if area of fracture changing state is relatively small
        msh = getMesh(obj,MortarSide.slave);
        areaChanged = sum(msh.surfaceArea(hasChangedElem));
        totArea = sum(msh.surfaceArea);
        if areaChanged/totArea < obj.activeSet.tol.areaTol
          obj.activeSet.curr = oldActiveSet;
          hasChanged = false;
          gresLog().log(1,['Active set update suppressed due to small fracture change:' ...
            ' areaChange/areaTot = %3.2e \n'],areaChanged/totArea);
        end
      end
    end



    function advanceState(obj)

      obj.stateOld = obj.state;
      obj.activeSet.prev = obj.activeSet.curr;

      % reset the counter for changed states
      obj.activeSet.stateChange(:) = 0;

      if obj.resetActiveSet
        % reset everything to stick for next time step
        obj.activeSet.curr(:) = ContactMode.stick;
      end
    end



    function goBackState(obj)
      % reset state to beginning of time step
      obj.state = obj.stateOld;

      if obj.resetActiveSet
        obj.activeSet.curr(:) = ContactMode.stick;
      else
        obj.activeSet.curr = obj.activeSet.prev;
      end

      obj.activeSet.stateChange(:) = 0;
    end



    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      outTraction = fac*obj.state.traction + ...
        (1-fac)*obj.stateOld.traction;
      outNormalGap = fac*obj.state.normalGap + ...
        (1-fac)*obj.stateOld.normalGap;
      outSlip = fac*obj.state.slip + ...
        (1-fac)*obj.stateOld.slip;
      outSliding = fac*obj.state.tangentialGap + ...
        (1-fac)*obj.stateOld.tangentialGap;

      tT = [outTraction(2:3:end),outTraction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      fractureState = ContactMode.integer(obj.activeSet.curr);

      pointStr = [];

      entries = {
        'normal_gap',              outNormalGap
        'slip_norm',               outSlip
        'normal_stress',           outTraction(1:3:end)
        'tangential_traction_1',   outTraction(2:3:end)
        'tangential_traction_2',   outTraction(3:3:end)
        'tangential_traction_norm',norm_tT
        'sliding_norm',            outSliding
        'fracture_state',          fractureState
        };

      surfaceStr = cell2struct(entries, {'name','data'}, 2);
    end

    function writeMatFile(obj,fac,tID)

      outTraction = fac*obj.state.traction + ...
        (1-fac)*obj.stateOld.traction;
      outNormalGap = fac*obj.state.normalGap + ...
        (1-fac)*obj.stateOld.normalGap;
      outSlip = fac*obj.state.slip + ...
        (1-fac)*obj.stateOld.slip;
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

      umOld = obj.domains(1).stateOld.data.displacements;
      usOld = obj.domains(2).stateOld.data.displacements;

      % stabilization contribution to the gap
      [H,~] = getStabilizationMatrixAndRhs(obj);
      stabGap = H*(obj.state.traction - obj.state.iniTraction);

      % recover variationally consistent gap (in local coordinates)
      currGap = obj.D*us + obj.M*um;
      prevGap = obj.D*usOld + obj.M*umOld;

      areaSlave = repelem(obj.interfMesh.msh(2).surfaceArea,3);

      gap = (currGap - stabGap)./areaSlave;
      slipIncrement = (currGap - prevGap - stabGap)./areaSlave;
      slipIncrement(1:3:end) = [];

      obj.state.normalGap = -gap(1:3:end);
      obj.state.tangentialGap = sqrt(gap(2:3:end).^2 + ...
        gap(3:3:end).^2);

      % norm of the slip
      obj.state.slip = sqrt(slipIncrement(1:2:end).^2 + ...
        slipIncrement(2:2:end).^2);

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
      um = obj.domains(1).state.data.displacements;
      us = obj.domains(2).state.data.displacements;
      um_old = obj.domains(1).stateOld.data.displacements;
      us_old = obj.domains(2).stateOld.data.displacements;

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


      for iPair = 1:obj.quadrature.numbInterfacePairs

        is = obj.quadrature.interfacePairs(iPair,1);
        im = obj.quadrature.interfacePairs(iPair,2);

        contactState = obj.activeSet.curr(is);

        % retrieve mortar integration data
        xiMaster = obj.quadrature.getMasterGPCoords(iPair);
        xiSlave = obj.quadrature.getSlaveGPCoords(iPair);
        dJw = obj.quadrature.getIntegrationWeights(iPair);

        % define slave related quantities
        nodeSlave = obj.interfMesh.local2glob{2}(obj.interfMesh.msh(2).surfaces(is,:));
        usDof = dofSlave.getLocalDoF(fldS,nodeSlave);
        tDof = getMultiplierDoF(obj,is);
        trac = obj.state.traction(tDof);
        % equilibrium equation and stabilization work with delta traction
        dTrac = trac - obj.state.iniTraction(tDof);

        % displacement increment to check if this is the first newton
        % iteration
        usDiff = norm(us(usDof) - us_old(usDof));


        nodeMaster = obj.interfMesh.local2glob{1}(obj.interfMesh.msh(1).surfaces(im,:));
        umDof = dofMaster.getLocalDoF(fldM,nodeMaster);

        [Nslave,Nmaster,Nmult] = ...
          getMortarBasisFunctions(obj.quadrature,im,is,xiMaster,xiSlave);

        % reshape basis function matrices to match number of components
        [Ns,Nm,Nmult] = reshapeBasisFunctions(3,Nslave,Nmaster,Nmult);

        % rotation of multipliers
        R = getRotationMatrix(obj.interfMesh,is);
        NmultR = pagemtimes(Nmult,R);

        % Reduce the dimension of multiplier basis functions exploiting
        % the local definition of degrees of freedom

        % the normal component of the multipliers basis
        Nmult_n = -Nmult(1,1,:);

        % tangential component of multiplier basis
        Nmult_t = NmultR(:,2:3,:);

        % get normal at the gauss points (for warped facets...?)
        normalNodes = obj.interfMesh.getNodalNormal(is);
        normal = pagemtimes(Ns,normalNodes);

        % operator selecting only tangential components of the
        % displacements
        T = eye(3) - pagemtimes(normal,'none',normal,'transpose');
        % normal and tangential component of displacement basis functions
        Nm_n = pagemtimes(normal,'transpose',Nm,'none');
        Nm_t = pagemtimes(T,Nm);

        Ns_n = pagemtimes(normal,'transpose',Ns,'none');
        Ns_t = pagemtimes(T,Ns);


        % tangential gap increment (quasi-static coulomb low)
        gt_curr = pagemtimes(Ns_t,us(usDof)) ...
          - pagemtimes(Nm_t,um(umDof));
        gt_old =  pagemtimes(Ns_t,us_old(usDof)) ...
          - pagemtimes(Nm_t,um_old(umDof));
        dgt = gt_curr - gt_old;


        % normal gap (scalar)
        g_n = pagemtimes(Ns_n,us(usDof)) ...
          - pagemtimes(Nm_n,um(umDof));


        % A_us
        Aun_m =  MortarQuadrature.integrate(f1,Nm_n,Nmult_n,dJw);
        Aun_s =  MortarQuadrature.integrate(f1,Ns_n,Nmult_n,dJw);
        asbMu.localAssembly(umDof,tDof(1),-Aun_m);
        asbDu.localAssembly(usDof,tDof(1),Aun_s);

        Aut_m =  MortarQuadrature.integrate(f1,Nm_t,Nmult_t,dJw);
        Aut_s =  MortarQuadrature.integrate(f1,Ns_t,Nmult_t,dJw);
        asbMu.localAssembly(umDof,tDof(2:3),-Aut_m);
        asbDu.localAssembly(usDof,tDof(2:3),Aut_s);

        % rhs (jump(eta),t)
        rhsUm(umDof) = rhsUm(umDof) - Aun_m*dTrac(1);
        rhsUs(usDof) = rhsUs(usDof) + Aun_s*dTrac(1);
        rhsUm(umDof) = rhsUm(umDof) - Aut_m*dTrac(2:3);
        rhsUs(usDof) = rhsUs(usDof) + Aut_s*dTrac(2:3);

        % assemble jacobian and rhs of traction balance equations

        % STICK MODE
        if contactState == ContactMode.stick

          asbMt.localAssembly(tDof(1),umDof,-Aun_m');
          asbDt.localAssembly(tDof(1),usDof,Aun_s');
          asbMt.localAssembly(tDof(2:3),umDof,-Aut_m');
          asbDt.localAssembly(tDof(2:3),usDof,Aut_s');

          % (mu_n,g_n)
          rhsT(tDof(1)) = rhsT(tDof(1)) ...
            + MortarQuadrature.integrate(f1,Nmult_n,g_n,dJw);
          % (mu_t,delta_g_t)
          rhsT(tDof(2:3)) = rhsT(tDof(2:3))  ...
            + MortarQuadrature.integrate(f1,Nmult_t,dgt,dJw);
        end

        % SLIP MODE
        if contactState == ContactMode.slip || contactState == ContactMode.newSlip

          slidingTol = obj.activeSet.tol.sliding;

          % A_nu
          Anu_m = MortarQuadrature.integrate(f1, Nmult_n, Nm_n, dJw);
          Anu_s = MortarQuadrature.integrate(f1, Nmult_n, Ns_n, dJw);
          asbMt.localAssembly(tDof(1),umDof,-Anu_m);
          asbDt.localAssembly(tDof(1),usDof,Anu_s);

          % A_tu (non linear term)
          if obj.state.slip(is) > slidingTol && usDiff > 100*eps
            % compute only on slip terms with sliding large enough
            dtdgt = computeDerTracGap(obj,trac(1),dgt);
            Atu_m = MortarQuadrature.integrate(f2, Nmult_t,dtdgt,Nm_t,dJw);
            Atu_s = MortarQuadrature.integrate(f2, Nmult_t,dtdgt,Ns_t,dJw);
            asbMt.localAssembly(tDof(2:3),umDof,Atu_m);
            asbDt.localAssembly(tDof(2:3),usDof,-Atu_s);

            % A_tn (non linear term)
            dtdtn = computeDerTracTn(obj,is,dgt);
            Atn = MortarQuadrature.integrate(f1,Nmult_t,dtdtn,dJw);
            asbQ.localAssembly(tDof(2:3),tDof(1),-Atn);
          end

          % A_tt
          Att = MortarQuadrature.integrate(f1,Nmult_t,Nmult_t,dJw);
          asbQ.localAssembly(tDof(2:3),tDof(2:3),Att);


          % rhs (mu_n,g_n)
          rhsT(tDof(1)) = rhsT(tDof(1)) + ...
            MortarQuadrature.integrate(f1,Nmult_n,g_n,dJw);

          % rhs -(mu_t,t*_T) (non linear term)
          tT_lim = computeLimitTraction(obj,is,dgt,trac,usDiff);

          %           % check deviation between limit traction and available traction
          %           dev = (R*trac)'*mean(tT_lim,3);
          %           if dev < 0
          %             %trust the limit traction
          %             trac(2:3) = -trac(2:3);
          %             obj.state.traction([3*is-1 3*is]) = -obj.state.traction([3*is-1 3*is]);
          %           end

          % rhs (mu_t,tT) tT and tT_lim in local coordinates (to be consistent with dof definition)
          rhsT(tDof(2:3)) = rhsT(tDof(2:3)) + Att*trac(2:3);
          % reconvert in local coordinates
          %             tT_lim = pagemtimes(R_el,'ctranspose',tT_lim,'none');
          rhsT(tDof(2:3)) = rhsT(tDof(2:3)) - ...
            MortarQuadrature.integrate(f1,Nmult_t,tT_lim,dJw);



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
          rhsT(tDof) = rhsT(tDof) + ...
            MortarQuadrature.integrate(f1,Nmult,trac,dJw);
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

      H = obj.stabilizationMat;

      elOpen = find(obj.activeSet.curr == ContactMode.open);
      elSlip = [find(obj.activeSet.curr == ContactMode.slip);...
        find(obj.activeSet.curr == ContactMode.newSlip)];

      dofOpen = DoFManager.dofExpand(elOpen,3);
      dofSlip = [3*elSlip-1; 3*elSlip];

      % remove rows and columns of dofs not requiring stabilization
      H([dofOpen;dofSlip],:) = 0;
      H(:,[dofOpen;dofSlip]) = 0;
      rhsH = -H*(obj.state.traction-obj.state.iniTraction);
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

    function dtdgt = computeDerTracGap(obj,sigma_n,dgt)
      % gt = obj.g_T(get_dof(nodeId));

      sz = size(dgt);
      dtdgt = zeros(3,3,sz(3));

      nG = sz(3);

      tauLim = obj.cohesion - tan(deg2rad(obj.phi))*sigma_n;

      for i = 1:nG  % gp loop
        g = dgt(:,:,i);
        dtdgt(:,:,i) = tauLim*((eye(3)*norm(g)^2 - g*g')/(norm(g))^3);
      end

    end

    function dtdtn = computeDerTracTn(obj,is,dgt)

      sz = size(dgt);
      dtdtn = zeros(sz);
      tanPhi = tan(deg2rad(obj.phi));

      if numel(sz)<4
        sz = [sz,1];
      end

      %if obj.slip.curr(is) > obj.activeSet.tol.sliding
      % use available gap to properly compute traction
      for i = 1:sz(3)
        g = dgt(:,:,i);
        dtdtn(:,:,i) = -tanPhi*(g/norm(g));
      end
      % else
      %         Rloc = obj.getRotationMatrix(is);
      %         t = Rloc*obj.state.traction(getMultiplierDoF(obj,is));
      %         % convert to global reference system
      %         dtdtn = repmat(-tanPhi*(t/norm(t)),1,1,sz(3),sz(4));
      % end
    end


    function tracLim = computeLimitTraction(obj,is,dgt,t,duNorm)

      % return the limit traction vector in the global frame
      sz = size(dgt);
      tracLim = zeros(sz);
      t_N = t(1);
      tauLim = obj.cohesion - tan(deg2rad(obj.phi))*t_N;


      if obj.state.slip(is) > obj.activeSet.tol.sliding && duNorm > 100*eps
        for i = 1:sz(3)
          g = dgt(:,:,i);
          tracLim(:,:,i) = tauLim*(g/norm(g));
        end
      else
        % compute tangential traction from traction
        R = obj.interfMesh.getRotationMatrix(is);
        Rt = R(:,2:3);
        t = Rt*t(2:3);
        tracLim =  repmat(tauLim*(t/norm(t)),1,1,sz(3));
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
        logical(getXMLData(input,0,"resetActiveSet"));
      contactSolver.activeSet.forceStickBoundary = ...
        logical(getXMLData(input,1,"forceStickBoundary"));

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
      contactSolver.activeSet.tol.maxStateChange = getXMLData(input,1e-4,"maxActiveSetChange");

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


