classdef SolidMechanicsContact < MeshTying

  % solid mechanics solver using piece-wise constant multipliers

  properties
    phi               % friction angle in deg
    cohesion          % cohesion
    contactHelper
    activeSet
  end


  methods
    function obj = SolidMechanicsContact(id,inputStruct,domains)

      obj@MeshTying(id,inputStruct,domains);

    end

    function registerInterface(obj,varargin)

      input = varargin{1};

      obj.stabilizationScale = getXMLData(input,1.0,"stabilizationScale");

      nDofsInterface = getNumbDoF(obj);

      obj.state.traction = zeros(nDofsInterface,1);
      obj.state.iniTraction = obj.state.traction;
      obj.state.slip = zeros(round(1/3*nDofsInterface),1);
      obj.state.normalGap = zeros(round(1/3*nDofsInterface),1);
      obj.state.tangentialGap = zeros(round(1/3*nDofsInterface),1);

      obj.stateOld = obj.state;

      initializeActiveSet(obj,input);

      computeConstraintMatrices(obj);

    end

    computeContactMatricesAndRhs(obj)

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

      computeContactMatricesAndRhs(obj);

      computeStabilizationMatrix(obj)

      % get stabilization matrix depending on the current active set
      [H,rhsStab] = getStabilizationMatrixAndRhs(obj);

      obj.Jconstraint = obj.Jconstraint - H;
      obj.rhsConstraint = obj.rhsConstraint + rhsStab;

      if gresLog().getVerbosity > 1
        % print rhs terms for each fracture state
        N = obj.nMult;
        dof_stick = DoFManager(obj.activeSet.curr == ContactMode.stick,3);
        dof_slip = [dofId(find(isSlip(obj.contactHelper,N')),3); dofId(find(isNewSlip(obj.contactHelper,N')),3)];
        dof_open = dofId(find(isOpen(obj.contactHelper,N')),3);
        fprintf('Rhs norm for stabilization: %4.3e \n', norm(rhsStab));
        fprintf('Rhs norm for stick dofs: %4.3e \n', norm(obj.rhsConstraint(dof_stick)))
        fprintf('Rhs norm for slip dofs: %4.3e \n', norm(obj.rhsConstraint(dof_slip)))
        fprintf('Rhs norm for open dofs: %4.3e \n', norm(obj.rhsConstraint(dof_open)))
      else
        return
      end

    end



    function goOnState(obj)
      % advance state to new time step
      obj.stateOld.traction = obj.state.traction;
      obj.normalGap.prev = obj.normalGap.curr;
      obj.tangentialGap.prev = obj.tangentialGap.curr;
      obj.slip.prev = obj.slip.curr;
      obj.contactHelper.activeSet.prev = obj.contactHelper.activeSet.curr;
      obj.contactHelper.activeSet.stateChange(:) = 0;
      if obj.contactHelper.resetActiveSet
        obj.contactHelper.activeSet.curr(:) = 1;   % reset everything to stick
      end
    end

    function goBackState(obj)
      % reset state to beginning of time step
      obj.state.traction = obj.stateOld.traction;
      obj.normalGap.curr = obj.normalGap.prev;
      obj.tangentialGap.curr = obj.tangentialGap.prev;
      obj.slip.curr = obj.slip.prev;
      if obj.contactHelper.resetActiveSet
        obj.contactHelper.activeSet.curr(:) = 1;
      else
        obj.contactHelper.activeSet.curr = obj.contactHelper.activeSet.prev;
      end
      obj.contactHelper.activeSet.stateChange(:) = 0;
    end


    function applyBCmaster(obj,bc,t)
      ph = obj.domains(1).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.domains(1).getSolver(ph),bc,t);
      obj.rhsUm(bcEnts) = 0;
      obj.Jmu(bcEnts,:) = 0;
      obj.Jmt(:,bcEnts) = 0;
    end


    function applyBCslave(obj,bc,t)
      ph = obj.domains(2).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.domains(2).getSolver(ph),bc,t);
      obj.rhsUs(bcEnts) = 0;
      obj.Jsu(bcEnts,:) = 0;
      obj.Jst(:,bcEnts) = 0;
    end


    function [cellStr,pointStr] = buildPrintStruct(obj,fac)

      nCellData = 8;

      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      if nargin ==1
        outTraction = obj.stateOld.traction;
        outNormalGap = obj.normalGap.prev;
        outSlip = obj.slip.prev;
        outSliding = obj.tangentialGap.prev;
      elseif nargin == 2
        outTraction = fac*obj.state.traction + ...
          (1-fac)*obj.stateOld.traction;
        outNormalGap = fac*obj.normalGap.curr + ...
          (1-fac)*obj.normalGap.prev;
        outSlip = fac*obj.slip.curr + ...
          (1-fac)*obj.slip.prev;
        outSliding = fac*obj.tangentialGap.curr + ...
          (1-fac)*obj.tangentialGap.prev;
      else
        error('Invalid number of input arguments for function buildPrintStruct');
      end

      % compute tangential norms
      tT = [outTraction(2:3:end),outTraction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      %       gapT = [outNormalGap(2:3:end),outNormalGap(3:3:end)];
      %       norm_gapT = sqrt(gapT(:,1).^2 + gapT(:,2).^2);

      fractureState = obj.contactHelper.activeSet.curr;
      %       fractureState(fractureState==3) = 2;
      %       fractureState(fractureState==4) = 3;

      cellStr(1).name = 'normal_gap';
      cellStr(1).data = outNormalGap;
      cellStr(2).name = 'slip_norm';
      cellStr(2).data = outSlip;
      cellStr(3).name = 'normal_stress';
      cellStr(3).data = outTraction(1:3:end);
      cellStr(4).name = 'tangential_traction_1';
      cellStr(4).data = outTraction(2:3:end);
      cellStr(5).name = 'tangential_traction_2';
      cellStr(5).data = outTraction(3:3:end);
      cellStr(6).name = 'tangential_traction_norm';
      cellStr(6).data = norm_tT;
      % print the fracture state
      cellStr(7).name = 'fracture_state';
      cellStr(7).data = fractureState; % [1 stick, 2 slip, 3 open]
      cellStr(8).name = 'sliding_norm';
      cellStr(8).data = outSliding; % [1 stick, 2 slip, 3 open]

      pointStr = [];        % when using P0 multipliers
    end



  end

  methods (Access = protected)

    function initializeActiveSet(obj,input)

      N = getMesh(obj,MortarSide.slave).nSurfaces;
      obj.activeSet.curr = repmat(ContactMode.stick,N,1);
      obj.activeSet.prev = obj.activeSet.curr;
      % count how many times an element has changed state during iteration
      obj.activeSet.stateChange = zeros(N,1);


      if isfield(input,"Tolerances")
        input = input.Tolerances;
      end

      obj.activeSet.tol.sliding = getXMLData(input,1e-4,"sliding");
      obj.activeSet.tol.normalGap = getXMLData(input,1e-6,"normalGap");
      obj.activeSet.tol.normalTrac = getXMLData(input,1e-3,"normalTraction");
      obj.activeSet.tol.slidingCheck = getXMLData(input,1e-4,"tangentialViolation");
      obj.activeSet.tol.minLimitTraction = getXMLData(input,1e-4,"minLimitTraction");
      obj.activeSet.tol.areaTol = getXMLData(input,1e-2,"areaChange");
      obj.activeSet.tol.maxStateChange = getXMLData(input,1e-4,"maxActiveSetChange");

    end


    function computeGap(obj)
      % compute normal gap and tangential slip (local coordinates)

      %nS = obj.interfMesh.msh(2).nSurfaces;

      um = obj.domains(1).state.data.displacements;
      us = obj.domains(2).state.data.displacements;

      umOld = obj.domains(1).stateOld.data.displacements;
      usOld = obj.domains(2).stateOld.data.displacements;

      % stabilization contribution to the gap
      [H,~] = getStabilizationMatrixAndRhs(obj);
      stabGap = H*(obj.state.traction - obj.iniTraction);

      % recover variationally consistent gap (in local coordinates)
      currGap = obj.D*us + obj.M*um;
      prevGap = obj.D*usOld + obj.M*umOld;

      areaSlave = repelem(obj.interfMesh.msh(2).surfaceArea,3);

      gap = (currGap - stabGap)./areaSlave;
      slipIncrement = (currGap - prevGap - stabGap)./areaSlave;
      slipIncrement(1:3:end) = [];

      obj.normalGap.curr = -gap(1:3:end);
      obj.tangentialGap.curr = sqrt(gap(2:3:end).^2 + ...
        gap(3:3:end).^2);

      % norm of the slip
      obj.slip.curr = sqrt(slipIncrement(1:2:end).^2 + ...
        slipIncrement(2:2:end).^2);

      % rotate gap from global to local coordinates
      %       for i = 1:nS
      %         R  = getRotationMatrix(obj.contactHelper,i);
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






    function hasChanged = updateActiveSet(obj)

      oldActiveSet = obj.contactHelper.activeSet.curr;

      for is = 1:numel(obj.contactHelper.activeSet.curr)

        % force boundary element to stick state
        nodes = obj.interfMesh.msh(2).surfaces(is,:);
        dofs = dofId(obj.interfMesh.local2glob{2}(nodes),3);
        if obj.contactHelper.forceStickBoundary
          if any(ismember(dofs,obj.dirDofs))
            % boundary element
            continue
          end
        end

        id = dofId(is,3);

        t = obj.state.traction(id);

        if obj.domains(2).simparams.verbosity > 2
          % report active set changes
          fprintf('\n Element %i: traction: %1.4e %1.4e %1.4e \n',is,t(:))
        end

        if isOpen(obj.contactHelper,is)
          if obj.normalGap.curr(is) < -  obj.contactHelper.tol.normalGap
            setStick(obj.contactHelper,is);
          end
        elseif t(1) > obj.contactHelper.tol.normalTrac
          setOpen(obj.contactHelper,is);
        else % not open
          tau = norm(t(2:3));
          tauLimit = abs(obj.cohesion - tan(deg2rad(obj.phi))*t(1));
          if tauLimit < obj.contactHelper.tol.minLimitTraction
            tauLimit = 0;
          end
          % relax stick/sliding transition
          if isStick(obj.contactHelper,is) && tau >= tauLimit
            tau = tau*(1-obj.contactHelper.tol.slidingCheck);
          elseif ~isStick(obj.contactHelper,is)  && tau <=tauLimit
            tau = tau*(1+obj.contactHelper.tol.slidingCheck);
          end
          if tau>tauLimit
            if isStick(obj.contactHelper,is)
              setNewSlip(obj.contactHelper,is);
            else
              setSlip(obj.contactHelper,is);
            end
          else
            setStick(obj.contactHelper,is);
          end
        end
      end

      % check if active set changed
      diffState = obj.contactHelper.activeSet.curr - oldActiveSet;
      idNewSlipToSlip = all([oldActiveSet==2 diffState==1],2);
      diffState(idNewSlipToSlip) = 0;
      hasChangedElem = diffState~=0;
      obj.contactHelper.activeSet.stateChange(hasChangedElem) = ...
        obj.contactHelper.activeSet.stateChange(hasChangedElem) + 1;

      hasChanged = any(diffState);

      if obj.domains(2).simparams.verbosity > 1
        % report active set changes
        as = obj.contactHelper.activeSet.curr;
        da = obj.contactHelper.activeSet.curr - oldActiveSet;
        d = da(oldActiveSet==1);
        assert(~any(d==2));
        fprintf('%i elements from stick to new slip \n',sum(d==1));
        fprintf('%i elements from stick to open \n',sum(d==3));
        d= da(oldActiveSet==2);
        fprintf('%i elements from new slip to stick \n',sum(d==-1));
        fprintf('%i elements from new slip to slip \n',sum(d==1));
        fprintf('%i elements from new slip to open \n',sum(d==2));
        d = da(oldActiveSet==3);
        fprintf('%i elements from slip to stick \n',sum(d==-2));
        fprintf('%i elements from slip to open \n',sum(d==1));
        d = da(oldActiveSet==4);
        fprintf('%i elements from open to stick \n',sum(d==-3));

        fprintf('Stick dofs: %i    Slip dofs: %i    Open dofs: %i \n',...
          sum(as==1), sum(any([as==2,as==3],2)), sum(as==4));
      end

      if hasChanged
        % check if area of fracture changing state is small
        areaChanged = sum(obj.interfMesh.msh(2).surfaceArea(hasChangedElem));
        totArea = sum(obj.interfMesh.msh(2).surfaceArea);
        if areaChanged/totArea < obj.contactHelper.tol.areaTol
          obj.contactHelper.activeSet.curr = oldActiveSet;
          hasChanged = false;
          if  obj.domains(2).simparams.verbosity > 1
            fprintf(['Active set update suppressed due to small fracture change:' ...
              ' areaChange/areaTot = %3.2e \n'],areaChanged/totArea);
          end
          return
        end
      end
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

      ncomp = dofSlave.getDoFperEnt(obj.physic);

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

      nDofMaster = dofMaster.getNumbDoF(obj.physic);
      nDofSlave = dofSlave.getNumbDoF(obj.physic);
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

      %if obj.slip.curr(is) > obj.contactHelper.tol.sliding
      % use available gap to properly compute traction
      for i = 1:sz(3)
        g = dgt(:,:,i);
        dtdtn(:,:,i) = -tanPhi*(g/norm(g));
      end
      % else
      %         Rloc = obj.contactHelper.getRotationMatrix(is);
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
  end

end


