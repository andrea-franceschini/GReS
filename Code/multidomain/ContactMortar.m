classdef ContactMortar < Mortar

  properties
    normalGap        % normal gap
    tangentialGap    % norm of tangential gap 
    slip             % norm of tangential slip
    dispJump         % displacement jump in global coordinates
    Jmu               % master displacement block
    Jmt               % master multipliers block
    Jsu               % slave displacement block
    Jst               % slave multipliers block
    Jtt               % traction - traction block
    rhsUm             % master displacement rhs
    rhsUs             % slave displacement rhs
    rhsT              % traction
    traction          % current and previous local tractions
    iniTraction     
    phi               % friction angle in deg
    cohesion   
    stabMat              % stabilization matrix
    stabScale
    %totMult
  end

  properties 
    contact
  end
  
  methods
    function obj = ContactMortar(id,inputStruct,domains)

      obj@Mortar(inputStruct,domains);

      obj.id = id;

      obj.multiplierType = "P0";
      obj.physic = "Poromechanics";

      isFldMaster = isField(obj.dofm(1),obj.physic);
      isFldSlave = isField(obj.dofm(2),obj.physic);
      assert(isFldMaster,['%s not available for ' ...
        'master domain %i'],obj.physic,obj.idDomain(1));
      assert(isFldSlave,['%s not available for ' ...
        'slave domain %i'],obj.physic,obj.idDomain(2));

      obj.phi = inputStruct.Coulomb.phiAttribute;
      obj.cohesion = inputStruct.Coulomb.cohesionAttribute;

      if isfield(inputStruct.Stabilization,"scaleAttribute")
        obj.stabScale = inputStruct.Stabilization.scaleAttribute;
      else
        obj.stabScale = 1.0;
      end

      %computing mortar matrices D and M and updating list of slave entitities
      computeMortarMatrices(obj);

      %remove inactive multipliers from D and M
      id = find(~any(obj.D,2));
      obj.D(id,:) = [];
      obj.M(id,:) = [];

      % initialize Jacobian and rhs for the interface (now that active
      % multipliers are known)
      initializeInterface(obj)

      finalizeInterface(obj.mesh,obj.solvers);
      
      % use the updated mesh object for the output of the interface
      if ~isempty(obj.outStruct)
        obj.outStruct.VTK = VTKOutput(obj.mesh.msh(2),obj.outStruct.name);
      end

      obj.contact = ContactHelper(obj);
    end

    function computeContactMortarMatrices(obj)

      % Compute contact matrices and rhs

      % the syntax of the comments matrices is inspired by the appendix of
      % Franceschini, Andrea, et al. "Algebraically stabilized Lagrange
      % multiplier method for frictional contact mechanics with
      % hydraulically active fractures." Computer Methods in Applied
      % Mechanics and Engineering 368 (2020): 113161.

      % extract current and previous displacements
      um = obj.solvers(1).state.data.dispCurr;
      us = obj.solvers(2).state.data.dispCurr;
      um_old = obj.solvers(1).state.data.dispConv;
      us_old = obj.solvers(2).state.data.dispConv;

      % compute displacement increment
%       um_slip = um - um_old;
%       us_slip = us - us_old;

      ncomp = obj.dofm(2).getDoFperEnt(obj.physic);

      % define matrix assemblers
      [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj);

      fldM = obj.dofm(1).getFieldId(obj.physic);
      fldS = obj.dofm(2).getFieldId(obj.physic);

      % anonymous functions for local fem computations
      f1 = @(a,b) pagemtimes(a,'ctranspose',b,'none');
      f2 = @(a,b,c) pagemtimes(a,'transpose',pagemtimes(b,c),'none');


      for is = 1:obj.mesh.msh(2).nSurfaces

        masterElems = find(obj.mesh.elemConnectivity(:,is));
        if isempty(masterElems)
          continue
        end

        % compute normal of each node in the element
        normal_nodes = obj.contact.getNodalNormal(is);

        % define slave related quantities
        nodeSlave = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(is,:));
        usDof = obj.dofm(2).getLocalDoF(nodeSlave,fldS);
        tDof = getMultiplierDoF(obj,is);
        trac = obj.traction.curr(tDof);
        % equilibrium equation and stabilization work with delta traction
        dTrac = trac - obj.iniTraction(tDof);

        % displacement increment to check if this is the first newton
        % iteration
        usDiff = norm(us(usDof) - us_old(usDof));
        
        %R_el = getRotationMatrix(obj.contact,is);

        for im = masterElems'

          nodeMaster = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(im,:));
          umDof = obj.dofm(1).getLocalDoF(nodeMaster,fldM);

          % compute basis function matrices on quadrature points
          [Nslave,Nmaster,Nmult] = ...
            getMortarBasisFunctions(obj.quadrature,is,im);

          if isempty(Nmaster)
            % refine connectivity matrix
            %obj.mesh.elemConnectivity(im,is) = 0;
            continue
          end

          % reshape basis function matrices to match number of components
          [Ns,Nm,Nmult] = ...
            obj.reshapeBasisFunctions(ncomp,Nslave,Nmaster,Nmult);

          % interpolate the normal at each gauss points (geometrical
          % approximation)
          normal = pagemtimes(Ns,normal_nodes);

          % get rotation matrix at each gauss point
          R_nod = getRotationMatrix(obj.contact,normal);

          % rotated multiplier basis functions
          NmultR = pagemtimes(Nmult,R_nod);
  
          % Reduce the dimension of multiplier basis functions exploiting
          % the local definition of degrees of freedome 

          % the normal component of the multipliers basis
          Nmult_n = -Nmult(1,1,:,:);      

          % tangential component of multiplier basis 
          Nmult_t = NmultR(:,2:3,:,:);

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
          Aun_m =  obj.quadrature.integrate(f1,Nm_n,Nmult_n);
          Aun_s =  obj.quadrature.integrate(f1,Ns_n,Nmult_n);
          asbMu.localAssembly(-Aun_m,umDof,tDof(1));
          asbDu.localAssembly(Aun_s,usDof,tDof(1));

          Aut_m =  obj.quadrature.integrate(f1,Nm_t,Nmult_t);
          Aut_s =  obj.quadrature.integrate(f1,Ns_t,Nmult_t);
          asbMu.localAssembly(-Aut_m,umDof,tDof(2:3));
          asbDu.localAssembly(Aut_s,usDof,tDof(2:3));

          % rhs (jump(eta),t)
          obj.rhsUm(umDof) = obj.rhsUm(umDof) - Aun_m*dTrac(1);
          obj.rhsUs(usDof) = obj.rhsUs(usDof) + Aun_s*dTrac(1);
          obj.rhsUm(umDof) = obj.rhsUm(umDof) - Aut_m*dTrac(2:3);
          obj.rhsUs(usDof) = obj.rhsUs(usDof) + Aut_s*dTrac(2:3);

          % assemble jacobian and rhs of traction balance equations

          % STICK MODE
          if isStick(obj.contact,is)

            asbMt.localAssembly(-Aun_m',tDof(1),umDof);
            asbDt.localAssembly(Aun_s',tDof(1),usDof);
            asbMt.localAssembly(-Aut_m',tDof(2:3),umDof);
            asbDt.localAssembly(Aut_s',tDof(2:3),usDof);

            % (mu_n,g_n)
            obj.rhsT(tDof(1)) = obj.rhsT(tDof(1)) ...
              + obj.quadrature.integrate(f1,Nmult_n,g_n);
            % (mu_t,delta_g_t)
            obj.rhsT(tDof(2:3)) = obj.rhsT(tDof(2:3))  ...
              + obj.quadrature.integrate(f1,Nmult_t,dgt);
          end

          % SLIP MODE
          if isSlip(obj.contact,is) || isNewSlip(obj.contact,is)

            % A_nu
            Anu_m = obj.quadrature.integrate(f1, Nmult_n, Nm_n);
            Anu_s = obj.quadrature.integrate(f1, Nmult_n, Ns_n);
            asbMt.localAssembly(-Anu_m,tDof(1),umDof);
            asbDt.localAssembly(Anu_s,tDof(1),usDof);

            % A_tu (non linear term)
            if obj.slip.curr(is) > obj.contact.tol.sliding && usDiff > 100*eps
              % compute only on slip terms with sliding large enough
              dtdgt = computeDerTracGap(obj,trac(1),dgt);
              Atu_m = obj.quadrature.integrate(f2, Nmult_t,dtdgt,Nm_t);
              Atu_s = obj.quadrature.integrate(f2, Nmult_t,dtdgt,Ns_t);
              asbMt.localAssembly(Atu_m,tDof(2:3),umDof);
              asbDt.localAssembly(-Atu_s,tDof(2:3),usDof);

              % A_tn (non linear term)
              dtdtn = computeDerTracTn(obj,is,dgt);
              Atn = obj.quadrature.integrate(f1,Nmult_t,dtdtn);
              asbQ.localAssembly(-Atn,tDof(2:3),tDof(1));
            end

            % A_tt
            Att = obj.quadrature.integrate(f1,Nmult_t,Nmult_t);
            asbQ.localAssembly(Att,tDof(2:3),tDof(2:3));

         
            % rhs (mu_n,g_n)
            obj.rhsT(tDof(1)) = obj.rhsT(tDof(1)) + ...
              obj.quadrature.integrate(f1,Nmult_n,g_n);


            % rhs (mu_t,tT) tT and tT_lim in local coordinates (to be consistent with dof definition) 
            obj.rhsT(tDof(2:3)) = obj.rhsT(tDof(2:3)) + Att*trac(2:3);

            % rhs -(mu_t,t*_T) (non linear term)
            tT_lim = computeLimitTraction(obj,is,dgt,trac,usDiff);
            % reconvert in local coordinates
%             tT_lim = pagemtimes(R_el,'ctranspose',tT_lim,'none');
            obj.rhsT(tDof(2:3)) = obj.rhsT(tDof(2:3)) - ...
              obj.quadrature.integrate(f1,Nmult_t,tT_lim);

            if obj.solvers(2).simparams.verbosity > 2
              fprintf('\nelement %i- rhsT: %5.3e %5.3e \n',is,Att*trac(2:3))
              fprintf('element %i- rhsTlim: %5.3e %5.3e \n',is,obj.quadrature.integrate(f1,Nmult_t,tT_lim))
              fprintf('------------------------------------ \n')
            end

          end

          % OPEN MODE
          if isOpen(obj.contact,is)

            % A_oo
            Aoo = obj.quadrature.integrate(f1,Nmult,Nmult);
            asbQ.localAssembly(Aoo,tDof,tDof);

            % rhs (mu,t)
            obj.rhsT(tDof) = obj.rhsT(tDof) + ...
              obj.quadrature.integrate(f1,Nmult,trac);
          end

        end % end inner master elems loop

      end   % end outer slave elems loop

      % assemble matrices into jacobian blocks
      obj.Jmu = asbMu.sparseAssembly();
      obj.Jsu = asbDu.sparseAssembly();
      obj.Jmt = asbMt.sparseAssembly();
      obj.Jst = asbDt.sparseAssembly();
      obj.Jtt = asbQ.sparseAssembly();

    end
    
    function initializeInterface(obj)
      nDofMult = getNumbMultipliers(obj);
      obj.dispJump = struct('prev',[],'curr',[]);
      obj.traction = struct('prev',[],'curr',[]);
      obj.slip = struct('prev',[],'curr',[]);
      obj.normalGap = struct('prev',[],'curr',[]);
      obj.traction.curr = zeros(nDofMult,1);
      obj.traction.prev = obj.traction.curr;
      obj.iniTraction = obj.traction.curr;
      obj.dispJump.curr = zeros(nDofMult,1);
      obj.dispJump.prev = obj.dispJump.curr;
      obj.slip.curr = zeros(round(1/3*nDofMult),1);
      obj.slip.prev = obj.slip.curr;
      obj.normalGap.curr = zeros(round(1/3*nDofMult),1);
      obj.normalGap.prev = obj.normalGap.curr;
      obj.totMult = nDofMult;
    end

    function resetRhs(obj)
      nDofMaster = getNumDoF(obj.dofm(1),obj.physic);
      nDofSlave = getNumDoF(obj.dofm(2),obj.physic);
      nDofMult = getNumbMultipliers(obj);
      obj.rhsUm = zeros(nDofMaster,1);
      obj.rhsUs = zeros(nDofSlave,1);
      obj.rhsT = zeros(nDofMult,1);
    end


    function computeMat(obj,~)
      resetRhs(obj);
      computeContactMortarMatrices(obj);
      computeStabilizationMatrix(obj);    % traction-jump stabilization
    end


    function computeRhs(obj)
      % the contact rhs is assembled directly in
      % computeMortarContactMatrices()
      [~,rhsStab] = getStabilizationMatrixAndRhs(obj); 
      obj.rhsT = obj.rhsT + rhsStab;
      if obj.solvers(2).simparams.verbosity > 1
        % print rhs terms for each fracture state
        N = 1:numel(obj.contact.activeSet.curr);
        dof_stick = dofId(find(isStick(obj.contact,N')),3);
        dof_slip = [dofId(find(isSlip(obj.contact,N')),3); dofId(find(isNewSlip(obj.contact,N')),3)]; 
        dof_open = dofId(find(isOpen(obj.contact,N')),3);
        fprintf('Rhs norm for stabilization: %4.3e \n', norm(rhsStab));
        fprintf('Rhs norm for stick dofs: %4.3e \n', norm(obj.rhsT(dof_stick)))
        fprintf('Rhs norm for slip dofs: %4.3e \n', norm(obj.rhsT(dof_slip)))
        fprintf('Rhs norm for open dofs: %4.3e \n', norm(obj.rhsT(dof_open)))
      else
        return
      end
    end


    function varargout = getJacobian(obj,field,domId)
      % get jacobian blocks associated to specific field and specific
      % domain
      % if nargout = 2 -> get master/slave pair of jacobian blocks
      % if nargout = 1 -> return multiplier jacobian

      if ~strcmp(field,obj.physic)
        % input field is not a mortar field
        varargout = cell(1,nargout);
        return
      end

      switch nargout
        case 1
          [H,~] = getStabilizationMatrixAndRhs(obj); 
          varargout{1} = obj.Jtt - H;
          % 0.5 is needed because the multiplier matrix is retrieved twice
        case 2
          % assign master/slave mortar matrix
          isDom = ismember(obj.idDomain,domId);
          if all(isDom)
            varargout{1} = obj.Jmu + obj.Jsu;
            varargout{2} = obj.Jmt + obj.Jst;
          elseif isDom(1)
            varargout{1} = obj.Jmu;
            varargout{2} = obj.Jmt;
          elseif isDom(2)
            varargout{1} = obj.Jsu;
            varargout{2} = obj.Jst;
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end


    function rhs = getRhs(obj,fld,varargin)
      % return rhs block associated to master/slave/multiplier field

      if ~strcmp(fld,obj.physic)
        rhs = [];
        return
      end

      switch nargin
        case 2
          % multiplier rhs block
          rhs = obj.rhsT;
        case 3
          domId = varargin{1};
          isDom = ismember(obj.idDomain,domId);
          if all(isDom)
            rhs = obj.rhsUm + obj.rhsUs;
          elseif isDom(1)
            rhs = obj.rhsUm;
          elseif isDom(2)
            rhs = obj.rhsUs;
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end


    function hasChanged = updateActiveSet(obj)

      oldActiveSet = obj.contact.activeSet.curr;

      for is = 1:numel(obj.contact.activeSet.curr)

        % force boundary element to stick state
        nodes = obj.mesh.msh(2).surfaces(is,:);
        dofs = dofId(obj.mesh.local2glob{2}(nodes),3);
        if obj.contact.forceStickBoundary
          if any(ismember(dofs,obj.dirDofs))
            % boundary element
            continue
          end
        end

        id = dofId(is,3);

        t = obj.traction.curr(id);

        if isOpen(obj.contact,is)
          if obj.dispJump(is) < -obj.contact.tol.normalGap
            setStick(obj.contact,is);
          end
        elseif t(1) > obj.contact.tol.normalTrac
          setOpen(obj.contact,is);
        else % not open
          tau = norm(t(2:3));
          tauLimit = abs(obj.cohesion - tan(deg2rad(obj.phi))*t(1));
          if tauLimit < obj.contact.tol.minLimitTraction
            tauLimit = 0;
          end
          % relax stick/sliding transition
          if isStick(obj.contact,is) && tau >= tauLimit
            tau = tau*(1-obj.contact.tol.slidingCheck);
          elseif ~isStick(obj.contact,is)  && tau <=tauLimit
            tau = tau*(1+obj.contact.tol.slidingCheck);
          end
          if tau>tauLimit
            if isStick(obj.contact,is)
              setNewSlip(obj.contact,is);
            else
              setSlip(obj.contact,is);
            end
          else
            setStick(obj.contact,is);
          end
        end
      end

      % check if active set changed
      diffState = obj.contact.activeSet.curr - oldActiveSet;
      idNewSlipToSlip = all([oldActiveSet==2 diffState==1],2);   
      diffState(idNewSlipToSlip) = 0;
      hasChangedElem = diffState~=0;
      hasChanged = any(diffState);
      if hasChanged
        % check if area of fracture changing state is small 
        areaChanged = sum(obj.mesh.msh(2).surfaceArea(hasChangedElem));
        totArea = sum(obj.mesh.msh(2).surfaceArea);
        if areaChanged/totArea < obj.contact.tol.areaTol
          hasChanged = false;
          if  obj.solvers(2).simparams.verbosity > 1
            fprintf(['Active set update suppressed due to small fracture change:' ...
              ' areaChange/areaTot = %3.2e \n'],areaChanged/totArea);
          end
          return
        end
      end

      if obj.solvers(2).simparams.verbosity > 1
        % report active set changes
        da = obj.contact.activeSet.curr - oldActiveSet;
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
      end
    end


    function updateState(obj,du)

      % traction update
      actMult = getMultiplierDoF(obj);
      n = numel(actMult);
      obj.traction.curr(actMult) = obj.traction.curr(actMult) + du(1:n);

      % update gap
      computeGap(obj);
    end

    function goOnState(obj)
      % advance state to new time step
      obj.traction.prev = obj.traction.curr;
      obj.dispJump.prev = obj.dispJump.curr;
      obj.slip.prev = obj.slip.curr;
    end

    function goBackState(obj)
      % reset state to beginning of time step
      obj.traction.curr = obj.traction.prev;
      obj.dispJump.curr = obj.dispJump.prev;
      obj.slip.curr = obj.slip.prev;
      obj.contact.activeSet.curr = obj.contact.activeSet.prev;
    end


    function applyBCmaster(obj,bc,t)
      ph = obj.solvers(1).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(ph),bc,t);
      obj.rhsUm(bcEnts) = 0;
      obj.Jmu(bcEnts,:) = 0;
      obj.Jmt(:,bcEnts) = 0;
    end


    function applyBCslave(obj,bc,t)
      ph = obj.solvers(2).bcs.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(ph),bc,t);
      bcEnts = removeSlaveBCdofs(obj,ph,bcEnts);
      obj.rhsUs(bcEnts) = 0;
      obj.Jsu(bcEnts,:) = 0;
      obj.Jst(:,bcEnts) = 0;
    end


    function [cellStr,pointStr] = buildPrintStruct(obj,fac)

      nCellData = 7;

      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      if nargin ==1
        outTraction = obj.traction.prev;
        outGap = obj.dispJump.prev;
        outSlip = obj.slip.prev;
      elseif nargin == 2
        outTraction = fac*obj.traction.curr + ...
          (1-fac)*obj.traction.prev;
        outGap = fac*obj.dispJump.curr + ...
          (1-fac)*obj.dispJump.prev;
        outSlip = fac*obj.slip.curr + ...
          (1-fac)*obj.slip.prev;
      else
        error('Invalid number of input arguments for function buildPrintStruct');
      end

      % compute tangential norms
      tT = [outTraction(2:3:end),outTraction(3:3:end)];
      norm_tT = sqrt(tT(:,1).^2 + tT(:,2).^2);

      gapT = [outGap(2:3:end),outGap(3:3:end)];
      norm_gapT = sqrt(gapT(:,1).^2 + gapT(:,2).^2);

      cellStr(1).name = 'normal_gap';
      cellStr(1).data = outGap(1:3:end);
      cellStr(2).name = 'tangential_gap_norm';
      cellStr(2).data = norm_gapT;
      cellStr(3).name = 'slip_norm';
      cellStr(3).data = outSlip;
      cellStr(4).name = 'normal_stress';
      cellStr(4).data = outTraction(1:3:end);
      cellStr(5).name = 'tangential_traction_1';
      cellStr(5).data = outTraction(2:3:end);
      cellStr(6).name = 'tangential_traction_2';
      cellStr(6).data = outTraction(3:3:end);
      cellStr(7).name = 'tangential_traction_norm';
      cellStr(7).data = norm_tT;

      pointStr = [];        % when using P0 multipliers
    end



  end

  methods (Access = protected)

    function [dofRow,dofCol,mat] = computeLoc(obj,kernel,dofRow,dofCol)

      if isnumeric(kernel)
        % local matrix is already provided
        mat = kernel;
      else
        % Compute local matrix using kernel as anonymous function
        mat = obj.quadrature.integrate(kernel);
      end

    end

    function computeGap(obj)
      % compute normal gap and tangential slip (local coordinates)

      nS = obj.mesh.msh(2).nSurfaces;

      um = obj.solvers(1).state.data.dispCurr;
      us = obj.solvers(2).state.data.dispCurr;

      umOld = obj.solvers(1).state.data.dispConv;
      usOld = obj.solvers(2).state.data.dispConv;

      % stabilization contribution to the gap
      [H,~] = getStabilizationMatrixAndRhs(obj);
      stabGap = H*(obj.traction.curr - obj.iniTraction);

      % recover variationally consistent gap
      currGap = obj.D*us + obj.M*um;
      prevGap = obj.D*usOld + obj.M*umOld;

      gap = (currGap - stabGap)./sum(obj.D,2);
      slipIncrement = (currGap - prevGap - stabGap)./sum(obj.D,2);

      % rotate gap from global to local coordinates
      for i = 1:nS
        R  = getRotationMatrix(obj.contact,i);
        obj.dispJump.curr(dofId(i,3)) = R'*gap(dofId(i,3));
        locSlip = R(:,2:3)'*slipIncrement(dofId(i,3));
        obj.slip.curr(i) = norm(locSlip);
      end


%       if obj.solvers(2).simparams.verbosity > 2
%         m = max(abs(stabGap./(obj.Jst*us)));
%         fprintf('Maximum gap deviation due to stabilization: %4.1f %% \n',100*m)
%       end

%       gapOld = (obj.D*us_old - obj.M*um_old - stabGap)./sum(obj.D,2);

    end

    function [H,rhsH] = getStabilizationMatrixAndRhs(obj)
      H = obj.stabMat;
      elList = (1:obj.mesh.msh(2).nSurfaces);
      % remove row/columns not requiring stabilization
      elOpen = find(isOpen(obj.contact,elList));
      dofOpen = dofId(elOpen,3);
      elSlip = [find(isSlip(obj.contact,elList));find(isNewSlip(obj.contact,elList))];
      dofSlip = [3*elSlip-1; 3*elSlip];
      H([dofOpen,dofSlip],:) = 0;
      H(:,[dofOpen,dofSlip]) = 0;
      rhsH = -H*(obj.traction.curr-obj.iniTraction);
    end





    function [asbMu,asbDu,asbMt,asbDt,asbQ] = defineAssemblers(obj)

      ncomp = obj.dofm(2).getDoFperEnt(obj.physic);

      % get number of index entries for sparse matrices

      % number of master nodes attacched to slave elements
      nNmaster = obj.mesh.msh(1).surfaceNumVerts'*obj.mesh.elemConnectivity;

      N1 = sum(nNmaster);
      N2 = sum(obj.mesh.elemConnectivity,1)*obj.mesh.msh(2).surfaceNumVerts;

      nmu = (ncomp^2)*N1;
      nsu = ncomp^2*N2;
      nmt = nmu;
      nst = nsu;
      nq = ncomp*N2;

      nDofMaster = obj.dofm(1).getNumDoF(obj.physic);
      nDofSlave = obj.dofm(2).getNumDoF(obj.physic);
      nDofMult = getNumbMultipliers(obj);

      % define matrix assemblers
      % anounymous function for local matrix assembly
      loc = @(kernel,dofR,dofC) computeLoc(obj,kernel,dofR,dofC);
      % initialize sparse matrix assemblers
      asbMu = assembler(nmu,loc,nDofMaster,nDofMult);
      asbDu = assembler(nsu,loc,nDofSlave,nDofMult);
      asbMt = assembler(nmt,loc,nDofMult,nDofMaster);
      asbDt = assembler(nst,loc,nDofMult,nDofSlave);
      asbQ = assembler(nq,loc,nDofMult,nDofMult);

    end

    function dtdgt = computeDerTracGap(obj,sigma_n,dgt)
      % gt = obj.g_T(get_dof(nodeId));

      sz = size(dgt);
      dtdgt = zeros(3,3,sz(3),sz(4));
      if numel(sz)<4
        sz = [sz,1];
      end

      nG = sz(3);

      tauLim = obj.cohesion - tan(deg2rad(obj.phi))*sigma_n;

      for i = 1:sz(4) % subtriangle loop
        for j = 1:nG  % gp loop
        g = dgt(:,:,j,i);
        dtdgt(:,:,j,i) = tauLim*((eye(3)*norm(g)^2 - g*g')/(norm(g))^3);
        end
      end

    end

    function dtdtn = computeDerTracTn(obj,is,dgt)

      sz = size(dgt);
      dtdtn = zeros(sz);
      tanPhi = tan(deg2rad(obj.phi));

      if numel(sz)<4
        sz = [sz,1];
      end

      %if obj.slip.curr(is) > obj.contact.tol.sliding
        % use available gap to properly compute traction
        for i = 1:sz(4)
          for j = 1:sz(3) % gp loop
            g = dgt(:,:,j,i);
            dtdtn(:,:,j,i) = -tanPhi*(g/norm(g));
          end
        end
     % else
%         Rloc = obj.contact.getRotationMatrix(is);
%         t = Rloc*obj.traction.curr(getMultiplierDoF(obj,is));
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

      if numel(sz)<4
        sz = [sz,1];
      end

      if obj.slip.curr(is) > obj.contact.tol.sliding && duNorm > 100*eps 
        for i = 1:sz(4) % sub triangle loop (for segment based)
          for j = 1:sz(3) % gp loop
            g = dgt(:,:,j,i);
            tracLim(:,:,j,i) = tauLim*(g/norm(g));
          end
        end
      else
        % compute tangential traction from traction
        R = obj.contact.getRotationMatrix(is);
        Rt = R(:,2:3);
        t = Rt*t(2:3);
        tracLim =  repmat(tauLim*(t/norm(t)),1,1,sz(3),sz(4));
      end
    end


    function computeStabilizationMatrix(obj)

      if ~isempty(obj.stabMat)
        % compute stabilization matrix only once for all edges
        % retrieve row-col needing stabilization at each time step
        return
      end

      % initialize matrix estimating number of entries
      % number of internal slave elements
      nes = sum(all(obj.mesh.e2f{2},2));
      nEntries = 2*36*nes; % each cell should contribute at least two times
      nmult = getNumbMultipliers(obj);
      localKernel = @(S,e1,e2) assembleLocalStabilization(obj,S,e1,e2);
      asbH = assembler(nEntries,localKernel,nmult,nmult);

      % get list of internal master edges
      inEdgeMaster = find(all(obj.mesh.e2f{1},2));

      for ieM = inEdgeMaster'
        % get master faces sharing internal edge ie
        fM = obj.mesh.e2f{1}(ieM,:);
        assert(numel(fM)==2,['Unexpected number of connected faces for' ...
          'master edge %i. Expected 2.'], ieM);

        % get slave faces sharing support with master faces
        fS = unique([find(obj.mesh.elemConnectivity(fM(1),:)),...
          find(obj.mesh.elemConnectivity(fM(2),:))]);

        if numel(fS) < 2
          continue
        end

        Am = mean(obj.mesh.msh(1).surfaceArea(fM));

        % get internal edges of slave faces

        eS = unique(obj.mesh.f2e{2}(fS,:));
        id = all(ismember(obj.mesh.e2f{2}(eS,:),fS),2);
        ieS = eS(id);

        % get active macroelement nodes
        nM = obj.mesh.e2n{1}(ieM,:);
        nS = unique(obj.mesh.e2n{2}(eS,:));

        % compute local schur complement approximation and retrieve
        % diagonal terms
        S = computeSchurLocal(obj,nM,nS,fS,obj.physic);

        % assemble stabilization matrix component
        for iesLoc = ieS'
          f = obj.mesh.e2f{2}(iesLoc,:);
          As = mean(obj.mesh.msh(2).surfaceArea(f));
          fLoc1 = dofId(find(fS==f(1)),3);
          fLoc2 = dofId(find(fS==f(2)),3);
          Sloc = 0.5*(Am/As)*(S(fLoc1,fLoc1)+S(fLoc2,fLoc2));
          asbH.localAssembly(Sloc,f(1),f(2));
        end
      end

      obj.stabMat = asbH.sparseAssembly();

      assert(norm(sum(obj.stabMat,2))<1e-8, 'Stabilization matrix is not locally conservative')
    end

    function S = computeSchurLocal(obj,nm,ns,fs,field)
      % compute approximate schur complement for local nonconforming
      % patch of element (in global reference)
      % input: nm/ns local master/slave node indices
      % fs: local slave faces indices

      % get slave and master dof to access jacobian
      fldM = getFieldId(obj.dofm(1),field);
      fldS = getFieldId(obj.dofm(2),field);
      dofS = obj.dofm(2).getLocalDoF(obj.mesh.local2glob{2}(ns),fldS);
      dofM = obj.dofm(1).getLocalDoF(obj.mesh.local2glob{1}(nm),fldM);


      % get local mortar matrices
      Dloc = obj.D(dofId(fs,3),dofS);
      Mloc = obj.M(dofId(fs,3),dofM);
      V = [Dloc, Mloc];              % minus sign!
      %V = Discretizer.expandMat(V,nc);

      % get local jacobian
      Km = getSolver(obj.solvers(1),field).J(dofM,dofM);
      Ks = getSolver(obj.solvers(2),field).J(dofS,dofS);
      Kloc = diag([1./diag(Ks);1./diag(Km)]);

      S = obj.stabScale*V*(Kloc*V');  % compute Schur complement
    end

    function [dofRow,dofCol,mat] = assembleLocalStabilization(obj,S,e1,e2)
      % assemble stabilization matrix S (in global coordinates) for
      % elements e1 and e2. 
      dof1 = dofId(e1,3);
      dof2 = dofId(e2,3);
      % get average rotation matrix
      n1 = getNormal(obj.contact,e1);
      n2 = getNormal(obj.contact,e2);
      if abs(n1'*n2 -1) < 1e4*eps
        avgR = ContactHelper.computeRot(n1);
      else
        A1 = obj.mesh.msh(2).surfaceArea(e1);
        A2 = obj.mesh.msh(2).surfaceArea(e2);
        nAvg = n1*A1 + n2*A2;
        nAvg = nAvg/norm(nAvg);
        avgR = ContactHelper.computeRot(nAvg);
      end

      % apply rotation matrix to S
      Srot = avgR'*S*avgR; 

      dofRow = [dof1;dof2];
      dofCol = [dof1;dof2];
      mat = [Srot,-Srot;-Srot,Srot];
      
    end


  end

end


