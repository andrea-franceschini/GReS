classdef EFEMaugmented < PhysicsSolver

  % solver for embedded tractions implementing the EFEM(0) formulation
  % Cusini et al (2021).

  properties
    activeSet = struct("curr",[],"prev",[])
    phi                     % friction angle in radians for each fracture
    cohesion                % the cohesion of each fracture
    fractureMesh            % a 2D mesh object with cut cell topology
    areaTol = 1e-6;         % minimum area of a fracture element
    bcTraction
    mechSolver              % handle to the Poromechanics solver
    fixActiveSet = false
    augmentation_t        % tangential augmentation coefficient cT used in t_T + cT*Delta g_T

  end

  properties (Access = private)
    fldMech
    fldFrac
    slipRecovery           % local recovery data for tangential slip DOFs
  end

  methods (Access = public)

    function obj = EFEMaugmented(domain)

      % call physicsSolver constructor
      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,varargin)

      obj.mechSolver = Poromechanics(obj.domain);
      obj.mechSolver.registerSolver(varargin{:});


      default = struct('Fracture',struct.empty,...
        'ActiveSet',missing,...
        'augmentation',1e3);

      params = readInput(default,varargin{:});
      obj.augmentation_t = params.augmentation;

      defineFractures(obj,params.Fracture);

      nCutCells = obj.fractureMesh.surfaces.num;

      dofm = obj.domain.dofm;

      % register nodal displacements on target regions
      dofm.registerVariable("fractureJump",entityField.cell,3,"nEntities",nCutCells);

      % store the id of the field in the degree of freedom manager
      flds = obj.getField;
      obj.fldMech = dofm.getVariableId(flds(1));
      obj.fldFrac = dofm.getVariableId(flds(2));

      % initialize the state object
      initState(obj);

      initializeActiveSet(obj,nCutCells,params.ActiveSet);

    end

    function assembleSystem(obj,dt)
      % compute the displacements matrices and rhs in the domain

      assembleSystemEFEM(obj,dt);

    end


    function initialize(obj)
      %
      obj.mechSolver.initialize()
      tIni = computeInitialTraction(obj);
      t = getState(obj,"traction");
      t = t + tIni(:);
      setState(obj,t,"traction");
      setStateOld(obj,t,"traction");
      setStateInit(obj,t,"traction");

    end

    function trac = computeInitialTraction(obj)
      % initialize traction for cell stress (average)
      avgStress = getState(obj.mechSolver,"avgStress");
      frac = obj.fractureMesh.surfaces;
      idx = [1;6;5;6;2;4;5;4;3];
      sigma = zeros(3);
      trac = zeros(3*frac.num,1);
      for i = 1:frac.num
        R = reshape(frac.rotationMatrices(i,:),3,3);
        cellId = frac.cutCells(i);
        sigma(:) = avgStress(cellId,idx);
        n = frac.normal(i,:);
        t = sigma*n';
        id = DoFManager.dofExpand(i,3);
        trac(id) = R'*t;
      end
    end


    function timeStepSetup(obj)

      % At the beginning of each time step, all non-open contact states are
      % reset to stick, except embedded faces where a traction boundary
      % condition is prescribed. Those faces are open from the start because
      % the EFEM residual must impose P*sigma = bcTraction by solving the
      % full u+w problem.
      hasTractionBC = obj.hasFractureTractionBC();
      isOpen = obj.activeSet.curr == ContactMode.open | hasTractionBC;
      obj.activeSet.curr(~isOpen) = ContactMode.stick;
      obj.activeSet.curr(hasTractionBC) = ContactMode.open;

    end


    function assembleSystemEFEM(obj,dt)
      % Assemble EFEM(0) system with state-dependent active fracture DOFs.
      %
      % Contact treatment implemented here:
      %   stick : solve only for continuous displacements u. The fracture
      %           jump w is anchored to zero and the traction is computed
      %           a posteriori from the projected 3D stress.
      %   slip  : solve globally for u and tangential jumps w_T1,w_T2,
      %           while the normal jump w_N is anchored to zero. Tractions
      %           are obtained from projected stresses and the Coulomb limit.
      %   open  : solve globally for u and all jump components w_N,w_T1,w_T2,
      %           imposing zero projected traction through the EFEM residual.

      % allocate
      dofm = obj.domain.dofm;
      frac = obj.fractureMesh.surfaces;
      cells = obj.grid.cells;

      subCells = dofm.getFieldCells(obj.fldMech);
      n = sum((obj.grid.nDim^2)*(obj.grid.cells.numVerts(subCells)).^2);
      n1 = sum((obj.grid.nDim^2)*(cells.numVerts(frac.cutCells)*frac.num));
      n2 = sum((obj.grid.nDim^2)*frac.num^2);
      nDofU = dofm.getNumbDoF(Poromechanics.getField());
      nDofW = dofm.getNumbDoF(obj.fldFrac);

      asbKuu = assembler(n,nDofU,nDofU);
      asbKuw = assembler(n1,nDofU,nDofW);
      asbKwu = assembler(n1,nDofW,nDofU);
      asbKww = assembler(n2,nDofW,nDofW);
      rhsU = zeros(nDofU,1);
      rhsW = zeros(nDofW,1);

      % local recovery containers for condensed tangential slip DOFs
      rec0 = struct('active',false,'uDof',[],'wNDof',[],'wTDof',[], ...
        'KTT',[],'KTx',[],'rT',[]);
      obj.slipRecovery = repmat(rec0,frac.num,1);

      % get state variables
      s = getState(obj);
      sOld = getStateOld(obj);
      time = s.time;
      iniStress = getStateInit(obj,'stress');
      iniTraction = getStateInit(obj,'traction');
      jump = s.fractureJump;
      jumpOld = sOld.fractureJump;
      du = s.displacements - sOld.displacements;
      dj = jump - jumpOld;

      coordinates = obj.grid.coordinates;

      % only hexa for now
      elem = Hexahedron(obj.grid,'gaussOrder',obj.mechSolver.getGaussOrder);
      nG = getNumbGaussPts(elem);

      cell2frac = zeros(cells.num,1);
      cell2frac(frac.cutCells) = 1:frac.num;

      topol = obj.grid.getCellNodes(subCells);

      for i = 1:numel(subCells)

        el = subCells(i);
        l = obj.mechSolver.cell2stress(el);

        f = cell2frac(el);
        isCellCut = f > 0;

        nodes = topol(i,:);
        uDof = dofm.getLocalDoF(obj.fldMech,nodes);
        coords = coordinates(nodes,:);

        % compute strain from continuous displacement field
        [gradN,dJw] = getDerBasisFAndDet(elem,coords);

        B = elem.getStrainMatrix(gradN);
        s.strain(l:l+nG-1,:) = reshape(pagemtimes(B,du(uDof)),6,nG)';

        if isCellCut

          % grab fracture DOFs and operators
          wDof = dofm.getLocalDoF(obj.fldFrac,f);
          Bw = computeCompatibilityMatrix(obj,frac,f,coords,gradN);
          E = computeEquilibriumOperator(obj,frac,f);

          state = obj.activeSet.curr(f);
          if obj.hasFractureTractionBC(f)
            state = ContactMode.open;
            obj.activeSet.curr(f) = ContactMode.open;
          end

          % Enhance strain only for components that are actually allowed in
          % the current contact mode. In stick the displacement is continuous.
          djEff = dj(wDof);
          if state == ContactMode.stick
            djEff(:) = 0;
            s.fractureJump(wDof) = 0;
          elseif state == ContactMode.slip || state == ContactMode.newSlip
            % Normal component is global. Tangential components are current
            % local/recovered values from the previous Newton iterate.
            % They are still included in the trial strain.
          elseif state == ContactMode.open
            % all components are global
          end

          enhancedStrain = reshape(pagemtimes(Bw,djEff),6,nG)';
          s.strain(l:l+nG-1,:) = s.strain(l:l+nG-1,:) + enhancedStrain;

        end

        % constitutive update
        [D, sigma, status] = obj.domain.materials.updateMaterial( ...
          cells.tag(el), ...
          sOld.stress(l:l+nG-1,:), ...
          s.strain(l:l+nG-1,:), ...
          dt, sOld.status(l:l+nG-1,:), el, time);

        % update stress map and gp counter
        s.status(l:l+nG-1,:) = status;
        s.stress(l:(l+nG-1),:) = sigma;

        % assemble internal forces for u
        dsigma = sigma - iniStress(l:l+nG-1,:);
        dsigma = reshape(dsigma',6,1,nG);
        sigma = reshape(sigma',6,1,nG); %#ok<NASGU>
        fTmp = pagemtimes(B,'ctranspose',dsigma,'none');
        fTmp = fTmp.*reshape(dJw,1,1,[]);
        fLoc = sum(fTmp,3);
        rhsU(uDof) = rhsU(uDof)+fLoc;

        KLoc = obj.mechSolver.computeKloc(B,D,B,dJw);
        asbKuu.localAssembly(uDof,uDof,KLoc);

        if isCellCut

          state = obj.activeSet.curr(f);
          if obj.hasFractureTractionBC(f)
            state = ContactMode.open;
            obj.activeSet.curr(f) = ContactMode.open;
          end
          fracId = frac.fracId(f);

          % full local EFEM operators, later reduced depending on state
          KuwLoc = Poromechanics.computeKloc(B,D,Bw,dJw);
          KwuLoc = Poromechanics.computeKloc(E,D,B,dJw);
          KwwLoc = Poromechanics.computeKloc(E,D,Bw,dJw);

          % projected stress residual on fracture
          fTmp = pagemtimes(E,'ctranspose',dsigma,'none');
          fTmp = fTmp.*reshape(dJw,1,1,[]);
          rSigma = sum(fTmp,3);

          tCurr = s.traction(wDof);
          tOld = sOld.traction(wDof);

          if state == ContactMode.stick

            % No fracture jump is solved in stick. Faces with prescribed
            % bcTraction have already been promoted to open, so this branch
            % only computes the reaction traction and anchors w.
            tracNew = iniTraction(wDof) + rSigma/frac.area(f) - obj.bcTraction(wDof);
            s.traction(wDof) = tracNew;

            asbKww.localAssembly(wDof,wDof,eye(3));
            rhsW(wDof) = rhsW(wDof) + s.fractureJump(wDof);

          elseif state == ContactMode.slip || state == ContactMode.newSlip

            % Slip: closed in the normal direction, free sliding in the
            % tangential directions. The condensed-tangential implementation
            % is not equivalent for the constant-sliding patch test because
            % it recovers w_T only locally.  Here w_T1,w_T2 are assembled as
            % global unknowns, while only w_N is anchored to zero.  This keeps
            % the no-penalty Coulomb law and restores global tangential
            % coupling along the fracture.
            projectedTrac = iniTraction(wDof) + rSigma/frac.area(f) - obj.bcTraction(wDof);
            [tracNew,slipDir,dDir,tauLim,tauSign] = updateTraction(obj,f,fracId,projectedTrac,tOld,tCurr, ...
              s.fractureJump(wDof),jumpOld(wDof));
            s.traction(wDof) = tracNew;

            rT = (tracNew - iniTraction(wDof))*frac.area(f);
            rBC = obj.bcTraction(wDof)*frac.area(f);
            rhsLoc = rSigma - rT - rBC;

            idN = 1;
            idT = [2 3];
            Afrac = frac.area(f);
            tanPhi = tan(obj.phi(fracId));

            % Tangential residual:
            %   R_T = rSigma_T - A*(t_T - tIni_T) - A*tBC_T,
            % with the augmented Coulomb direction
            %   t_T = tauLim*qhat,
            %   qhat = q/||q||,
            %   q = projectedTrac_T + cT*(g_T - g_T_old),
            %   tauLim = abs(c - tan(phi)*t_N),
            %   t_N = tIni_N + rSigma_N/A - tBC_N.
            % Since d(projectedTrac_T) = d(rSigma_T)/A,
            %   dq = d(rSigma_T)/A + cT*d(g_T).
            % Therefore
            %   dR_T = d(rSigma_T)
            %        + tan(phi)*sign(c - tan(phi)*t_N)*qhat*d(rSigma_N)
            %        - A*tauLim*dqhat,
            % with dqhat = Dq*dq.
            Iaug = eye(2) - tauLim*dDir;
            KTT = Iaug*KwwLoc(idT,idT) + tanPhi*tauSign*(slipDir*KwwLoc(idN,idT)) ...
                  - Afrac*tauLim*obj.augmentation_t*dDir;
            KTu = Iaug*KwuLoc(idT,:) + tanPhi*tauSign*(slipDir*KwuLoc(idN,:));
            KTN = Iaug*KwwLoc(idT,idN) + tanPhi*tauSign*(slipDir*KwwLoc(idN,idN));

            % Mechanical equation: depends on the active tangential jump.
            asbKuw.localAssembly(uDof,wDof(idT),KuwLoc(:,idT));

            % Tangential equilibrium equations.
            asbKwu.localAssembly(wDof(idT),uDof,KTu);
            asbKww.localAssembly(wDof(idT),wDof(idT),KTT);
            asbKww.localAssembly(wDof(idT),wDof(idN),KTN);
            rhsW(wDof(idT)) = rhsW(wDof(idT)) + rhsLoc(idT);

            % Normal contact closure: w_N = 0. The normal projected stress
            % still enters the Coulomb limit through the consistent tangent
            % above, but the normal residual itself is not assembled because
            % it is identically zero when t_N is defined from P_N*sigma.
            asbKww.localAssembly(wDof(idN),wDof(idN),1);
            rhsW(wDof(idN)) = rhsW(wDof(idN)) + s.fractureJump(wDof(idN));

          elseif state == ContactMode.open

            % Open: full displacement jump is global, and projected traction
            % is zero through updateTraction.
            projectedTrac = iniTraction(wDof) + rSigma/frac.area(f) - obj.bcTraction(wDof);
            tracNew = zeros(3,1);
            dtdg = zeros(3);
            s.traction(wDof) = tracNew;

            KwwLoc = KwwLoc - dtdg*frac.area(f);

            asbKuw.localAssembly(uDof,wDof,KuwLoc);
            asbKwu.localAssembly(wDof,uDof,KwuLoc);
            asbKww.localAssembly(wDof,wDof,KwwLoc);

            rT = (tracNew - iniTraction(wDof))*frac.area(f);
            rBC = obj.bcTraction(wDof)*frac.area(f);
            rhsLoc = rSigma - rT - rBC;
            rhsW(wDof) = rhsW(wDof) + rhsLoc;

          else
            error('Unsupported contact state in EFEM assembly')
          end

        end

      end

      % update modified state
      setState(obj,s);

      % populate rhs and jacobian
      obj.domain.J{obj.fldMech,obj.fldMech} = asbKuu.sparseAssembly;
      obj.domain.J{obj.fldMech,obj.fldFrac} = asbKuw.sparseAssembly;
      obj.domain.J{obj.fldFrac,obj.fldMech} = asbKwu.sparseAssembly;
      obj.domain.J{obj.fldFrac,obj.fldFrac} = asbKww.sparseAssembly;
      obj.domain.rhs{obj.fldMech} = rhsU;
      obj.domain.rhs{obj.fldFrac} = rhsW;

    end


    function initState(obj)
      % add embedded fracture fields to state structure
      nCutCells = obj.fractureMesh.surfaces.num;
      state = getState(obj);
      state.fractureJump = zeros(3*nCutCells,1);
      state.traction = zeros(3*nCutCells,1);

      % elastic and plastic slip are obtained subtracting current and
      % previous elastic and plastic jump
      setState(obj,state);
      obj.bcTraction = zeros(3*nCutCells,1);
      obj.activeSet.curr = repmat(ContactMode.stick,nCutCells,1);
      obj.activeSet.prev = obj.activeSet.curr;

    end


     function hasConfigurationChanged = updateConfiguration(obj)

      oldActiveSet = obj.activeSet.curr;

      traction = getState(obj,"traction");
      displacementJump = getState(obj,"fractureJump");

      f = obj.fractureMesh.surfaces;
      hasTractionBC = obj.hasFractureTractionBC();
      obj.activeSet.curr(hasTractionBC) = ContactMode.open;


      for i = 1:numel(obj.activeSet.curr)

        if hasTractionBC(i)
          obj.activeSet.curr(i) = ContactMode.open;
          continue
        end

        fracId = f.fracId(i);

        state = obj.activeSet.curr(i);

        id = DoFManager.dofExpand(i,3);

        t = traction(id);
        g_n = displacementJump(id(1));

        limitTraction = abs(obj.cohesion(fracId) - tan(obj.phi(fracId))*t(1));
        tangentialTraction = norm(t(2:3));
        yieldTol = 1e-10*max([1,limitTraction,tangentialTraction]);

        % No-penalty slip puts the traction exactly on the Coulomb surface.
        % Near equality, roundoff can place one element infinitesimally
        % inside the cone and make updateContactState classify it as stick.
        % Preserve sliding history on the yield surface, and allow stick to
        % enter slip when it is numerically on the surface.
        if tangentialTraction >= limitTraction - yieldTol
          if state == ContactMode.stick
            obj.activeSet.curr(i) = ContactMode.newSlip;
          elseif state == ContactMode.newSlip || state == ContactMode.slip
            obj.activeSet.curr(i) = ContactMode.slip;
          else
            obj.activeSet.curr(i) = updateContactState(state,t,...
                                    limitTraction, ...
                                    g_n,...
                                    obj.activeSet.tol);
          end
        else
          obj.activeSet.curr(i) = updateContactState(state,t,...
                                  limitTraction, ...
                                  g_n,...
                                  obj.activeSet.tol);
        end

      end

      % check if active set changed
      asNew = obj.activeSet.curr;
      asNew(hasTractionBC) = ContactMode.open;
      asOld = oldActiveSet;

      % do not upate state of element that exceeded the maximum number of
      % individual updates
      reset = obj.activeSet.stateChange >= ...
        obj.activeSet.tol.maxStateChange;
      reset(hasTractionBC) = false;

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
        areaChanged = sum(f.area(hasChangedElem));
        totArea = sum(f.area);
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

    function isReset = resetConfiguration(obj)

      hasTractionBC = obj.hasFractureTractionBC();
      toReset = obj.activeSet.curr(:) ~= ContactMode.open & ~hasTractionBC;
      obj.activeSet.curr(toReset) = ContactMode.stick;
      obj.activeSet.curr(hasTractionBC) = ContactMode.open;

      isReset = true;

    end


    function advanceState(obj)
      % Set converged state to current state after newton convergence

      % add elastic/plastic jump contribution to stored variables
      obj.mechSolver.advanceState();

      obj.activeSet.prev = obj.activeSet.curr;

    end


    function [tractionNew, slipDir, dDir, tauLim, tauSign] = updateTraction(obj,fEl,fracId,projectedTrac,tOld,tCurr,jumpNew,jumpOld) %#ok<INUSL>

      % projectedTrac is the traction obtained directly from the projected
      % 3D stress, in the local fracture frame.
      %
      % Local frame convention:
      %   component 1     : normal
      %   components 2, 3 : tangential directions

      tractionNew = zeros(3,1);
      slipDir = zeros(2,1);
      dDir = zeros(2);
      tauLim = 0;
      tauSign = 1;

      state = obj.activeSet.curr(fEl);

      if state == ContactMode.open
        % Open fracture: projected traction is zero.
        return
      end

      if state == ContactMode.stick
        % In stick, the jump is constrained to zero and traction is computed
        % a posteriori as the projection of the current 3D stress.
        tractionNew = projectedTrac;
        return
      end

      if state == ContactMode.slip || state == ContactMode.newSlip

        % Normal traction is not regularized: it comes from P_N*sigma.
        tractionNew(1) = projectedTrac(1);

        tauRaw = obj.cohesion(fracId) - tan(obj.phi(fracId))*tractionNew(1);
        tauLim = abs(tauRaw);
        if tauRaw < 0
          tauSign = -1;
        else
          tauSign = 1;
        end

        slip = jumpNew([2;3]) - jumpOld([2;3]);

        % Augmented sliding direction.  This is the only difference with the
        % penalty-free slip projection: the direction is not based only on
        % Delta g_T, but on q = t_T^proj + cT*Delta g_T.
        qAug = projectedTrac(2:3) + obj.augmentation_t*slip;
        qNorm = norm(qAug);
        isDirectionReliable = qNorm > obj.activeSet.tol.sliding;

        if isDirectionReliable
          slipDir = qAug/qNorm;
          dDir = (qNorm^2*eye(2) - qAug*qAug')/qNorm^3;
        else
          % If the augmented direction is too small, fall back to an
          % available tangential traction direction to avoid division by zero.
          if norm(tCurr(2:3)) > 0
            slipDir = tCurr(2:3)/norm(tCurr(2:3));
          elseif norm(tOld(2:3)) > 0
            slipDir = tOld(2:3)/norm(tOld(2:3));
          elseif norm(projectedTrac(2:3)) > 0
            slipDir = projectedTrac(2:3)/norm(projectedTrac(2:3));
          else
            slipDir = [1;0];
          end
          dDir = zeros(2);
        end

        tractionNew(2:3) = tauLim * slipDir;

      else
        error('Unsupported contact state in traction update')
      end

    end


    function order = getGaussOrder(obj)
      order = obj.mechSolver.getGaussOrder;
    end



    function goBackState(obj)

      % reset state to beginning of time step
      obj.mechSolver.goBackState();



      obj.activeSet.curr = obj.activeSet.prev;
      obj.activeSet.curr(obj.hasFractureTractionBC()) = ContactMode.open;

      if obj.activeSet.resetActiveSet
        resetConfiguration(obj);
      end

    end

    function updateState(obj,solution)

      % Update state structure with last solution increment.
      %
      % The DoFManager still owns three fractureJump components per cut cell.
      % Components removed from the global EFEM solve are anchored in the
      % assembled linear system. For slip, tangential components are updated
      % here from the local condensation recovery equation.
      obj.mechSolver.updateState(solution);

      dofm = obj.domain.dofm;
      ents = dofm.getActiveEntities(obj.fldFrac,1);
      stateCurr = obj.getState();

      if nargin > 1
        dw = solution(getDoF(dofm,obj.fldFrac));
        stateCurr.fractureJump(ents) = stateCurr.fractureJump(ents) + dw;

        % Enforce stick continuity explicitly.
        stick = find(obj.activeSet.curr == ContactMode.stick);
        if ~isempty(stick)
          dofStick = getLocalDoF(dofm,obj.fldFrac,stick);
          stateCurr.fractureJump(dofStick) = 0;
        end

        % Tangential slip DOFs are solved globally in slip. No local
        % recovery is needed here.

        setState(obj,stateCurr);
      end

    end


    function applyBC(obj,bcId,t)

      obj.mechSolver.applyBC(bcId,t);

    end

    function applyDirVal(obj,bcId,t)

      obj.mechSolver.applyDirVal(bcId,t);

    end


    function out = isLinear(obj)
      out = false;
    end

    function out = isSymmetric(obj)

      % if the problem is linear, then Poromechanics is symmetric
      out = false;

    end

    function writeSolution(obj,fac,tID)

      jumpOld = getStateOld(obj,"fractureJump");
      jumpCurr = getState(obj,"fractureJump");

      tractionOld = getStateOld(obj,"traction");
      tractionCurr = getState(obj,"traction");

      obj.domain.outstate.results(tID).fractureJump = jumpCurr*fac+jumpOld*(1-fac);
      obj.domain.outstate.results(tID).traction = tractionCurr*fac+tractionOld*(1-fac);

    end

    function [cellData,pointData] = writeVTK(obj,fac,time)

      [cellData,pointData] = obj.mechSolver.writeVTK(fac,time);


      cellDataFrac = repmat(struct('name', 1, 'data', 1), 1, 1);
      cellDataFrac(1).name = 'isCellFractured';
      isCutCell = false(obj.grid.cells.num,1);
      isCutCell(obj.fractureMesh.surfaces.cutCells) = true;
      cellDataFrac(1).data = double(isCutCell);

      cellData = OutState.mergeOutFields(cellData,cellDataFrac);


      % this method do not return outputs for the 3D mesh grid. instead, it
      % works on a separate 2D polygonal mesh
      s = obj.domain.state.interpolate(fac);
      jump = s.fractureJump;
      jumpDiff = s.fractureJump - obj.domain.getStateOld.fractureJump;
      trac = s.traction(:);

      nCellData = 3;
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      % Displacements
      cellStr(1).name = 'fractureJump';
      cellStr(1).data = [jump(1:3:end), jump(2:3:end), jump(3:3:end)];

      cellStr(2).name = 'traction';
      cellStr(2).data = [trac(1:3:end), trac(2:3:end), trac(3:3:end)];

      cellStr(3).name = 'fractureState';
      as = double(obj.activeSet.curr);
      cellStr(3).data = as;

      cellStr(4).name = 'jumpSlip';
      cellStr(4).data = [jumpDiff(2:3:end), jumpDiff(3:3:end)];

      % plot directly into the domain vtm block
      blk = obj.domain.vtmBlock;
      obj.domain.outstate.writeVTKfile(blk,'EmbeddedFractures',obj.fractureMesh,...,
        time,[],[],[],cellStr)

    end


  end




  methods (Access=private)

    function hasBC = hasFractureTractionBC(obj,fracId)
      % Return true for fracture-surface elements with a prescribed local
      % traction. A nonzero bcTraction means the element must be treated as
      % open from the active-set point of view, so the full u+w EFEM
      % residual imposes P*sigma = bcTraction.
      if isempty(obj.bcTraction)
        n = obj.fractureMesh.surfaces.num;
        hasAll = false(n,1);
      else
        bc = reshape(obj.bcTraction,3,[])';
        hasAll = any(abs(bc) > 0,2);
      end

      if nargin > 1
        hasBC = hasAll(fracId);
      else
        hasBC = hasAll;
      end
    end

    function defineFractures(obj,fractureStruct)

      % define the fracture geometrical informations

      coords = obj.grid.coordinates;

      nFractures = numel(fractureStruct);

      if nFractures == 0
        error("At least one fracture must be defined")
      end

      [edgesTopol,cellToEdges] = getEdgeTopology(obj);

      nE = size(edgesTopol,1);

      tol = 1e-8;

      obj.cohesion = zeros(nFractures,1);
      obj.phi = zeros(nFractures,1);

      obj.fractureMesh = Grid();
      fMesh = obj.fractureMesh;
      f = fMesh.surfaces;
      [f.tang1,f.tang2,f.cutCells,f.fracId] = deal([]);
      f.connectivity = ArrayOfArrays();

      % count number of vertices in each fracture
      nV = 0;

      for fId = 1:nFractures

        % read the fracture geometry
        d = double.empty;
        default = struct('normal',d,'origin',d,...
          'dimensions',d,'lengthVec',d,'widthVec',d,...
          'cohesion',d,'frictionAngle',d);

        frac = readInput(default,fractureStruct(fId));

        normal = frac.normal;
        origin = frac.origin;
        dims = frac.dimensions;
        lVec = frac.lengthVec;
        wVec = frac.widthVec;
        obj.cohesion(fId) = frac.cohesion;
        obj.phi(fId) = deg2rad(frac.frictionAngle);

        assert(all([abs(lVec*normal') < 1e-3,...
          abs(wVec*normal')<1e-3,...
          abs(lVec*wVec')<1e-3]),...
          "For now, normal, length and width direction must form " + ...
          "an orthonormal basis");

        % define plane corners
        lVec = lVec/norm(lVec);   wVec = wVec/norm(wVec);
        L = 0.5*dims(1)*lVec;
        W = 0.5*dims(2)*wVec;
        A = origin - L - W;
        B = origin + L - W;
        C = origin + L + W;
        D = origin - L + W;

        % mark node location w.r.t plane
        distVec = coords - origin;
        nVec = reshape(normal/norm(normal),[],1);
        R = mxComputeRotationMat(nVec);
        tVec1 = R(:,2);
        tVec2 = R(:,3);

        sign = distVec * nVec > 0;

        % check edges with nodes in different locations
        isEdgeCrossed = sum(sign(edgesTopol),2) == 1;

        % compute intersection points of crossed edges
        x0 = coords(edgesTopol(isEdgeCrossed,1),:);
        x1 = coords(edgesTopol(isEdgeCrossed,2),:);
        distEdge = x1-x0;
        distOrigin = origin - x0;
        t = (distOrigin * nVec)./(distEdge * nVec);
        xInt = x0 + t.*distEdge;

        intersections = zeros(nE,3);
        intersections(isEdgeCrossed,:) = xInt;

        % check if points lie in the finite plane
        abVec = B - A;
        adVec = D - A;
        vec = xInt - A;
        isInPlane1 = abs((xInt - origin) * nVec) < tol;
        dotAB = vec * abVec';
        dotAD = vec * adVec';
        isInPlane2 = dotAB > -tol & dotAB < norm(abVec)^2 + tol;
        isInPlane3 = dotAD > -tol & dotAD < norm(adVec)^2 + tol;

        isInPlane = all([isInPlane1,isInPlane2,isInPlane3],2);

        crossedEdgeId = find(isEdgeCrossed);
        tipEdges = crossedEdgeId(~isInPlane);

        edgeTag = zeros(length(isEdgeCrossed),1);
        edgeTag(isEdgeCrossed) = 1;
        edgeTag(tipEdges) = 2;

        % mark cut cells and update object properties
        m = reshape(edgeTag(cellToEdges),[],size(cellToEdges,2));
        isCutCell = any(m,2) & ~any(m==2,2);

        nC = sum(isCutCell);

        newCutCells = find(isCutCell);

        % preallocate number of nodes (maximum 6 per cell in hexa)
        surfs = zeros(nC*6,1);
        cutCellVertices = zeros(nC*6,3);
        cutNumVerts = zeros(nC,1);

        % loop over cut cell and compute geometry
        c = 0;

        for ic = 1:nC
          cellEdges = cellToEdges(newCutCells(ic),:);
          isEdgeCut = logical(m(newCutCells(ic),:));
          cutEdges = cellEdges(isEdgeCut);
          nVerts = numel(cutEdges);
          verts = intersections(cutEdges,:);
          [cutCellVertices(c+1:c+nVerts,:),perm] = orderPointsCCW(verts);
          cutNumVerts(ic) = nVerts;
          surfs(c+1:c+nVerts) = cutEdges(perm);
          c = c + nVerts;
        end

        cutCellVertices = cutCellVertices(1:sum(cutNumVerts),:);

        % compute geometry
        normals = repmat(nVec',nC,1);
        [cutA,cutC] = computePolygonGeometry(cutCellVertices,cutNumVerts,normals);

        f.cutCells    = [f.cutCells; newCutCells];
        f.normal      = [f.normal; normals];
        f.area        = [f.area; cutA];
        f.center      = [f.center; cutC];
        f.tang1       = [f.tang1; repmat(tVec1',nC,1)];
        f.tang2       = [f.tang2; repmat(tVec2',nC,1)];
        f.fracId  = [f.fracId; repmat(fId,nC,1)];

        % finalize the mesh for the current fracture
        surfs = surfs(1:c);

        % local numbering w.r.t. current fracture
        [usedEdges,~,locConn] = unique(surfs,'stable');

        % shift to global numbering in the accumulated fracture mesh
        surfs = nV + locConn;
        nV = nV + numel(usedEdges);

        f.connectivity = [f.connectivity; ArrayOfArrays(surfs,cutNumVerts)];
        f.numVerts = [f.numVerts; cutNumVerts];

        % coordinates consistent with local numbering above
        fMesh.coordinates = [fMesh.coordinates; intersections(usedEdges,:)];

      end

      % discard too small fractures
      id = f.area > obj.areaTol;
      f.num           = sum(id);
      f.cutCells      = f.cutCells(id);
      f.center        = f.center(id,:);
      f.area          = f.area(id);
      f.normal        = f.normal(id,:);
      f.tang1         = f.tang1(id,:);
      f.tang2         = f.tang2(id,:);
      f.numVerts      = f.numVerts(id);
      f.VTKType       = repmat(double(VTKType.Polygon),f.num,1);

      surfs = getRows(f.connectivity,find(id));

      [u,~,id2] = unique(getData(surfs));
      f.connectivity = ArrayOfArrays(id2(:),f.numVerts);
      fMesh.coordinates = fMesh.coordinates(u,:);

      R = zeros(f.num,9);
      for i = 1:f.num
        n = f.normal(i,:);
        Ri = mxComputeRotationMat(n);
        R(i,:) = Ri(:);
      end
      f.rotationMatrices = R;

      % finalize the grid
      fMesh.surfaces = f;
      %initializeGrid(fMesh);
      obj.fractureMesh = fMesh;

    end


    function [edges,c2e] = getEdgeTopology(obj)

      cells = obj.grid.cells;
      vtkId = unique(cells.VTKType);
      assert(isscalar(vtkId),"EFEM implemented only for mesh with uniform element shapes");

      if vtkId == VTKType.Hexa

        eLoc = [ ...
          1 2; 2 3; 3 4; 4 1; ...
          5 6; 6 7; 7 8; 8 5; ...
          1 5; 2 6; 3 7; 4 8 ];

        nCells = cells.num;
        nEdgesLoc = size(eLoc,1);

        topol = obj.grid.getCellNodes(1:nCells);

        % Node indices of all edges (stacked)
        e1 = topol(:, eLoc(:,1));   % (nCells x 12)
        e2 = topol(:, eLoc(:,2));   % (nCells x 12)

        allEdges = [e1(:), e2(:)];  % (12*nCells x 2)
        allEdges = sort(allEdges, 2);

        % discard duplicated edges
        [edges, ~, ic] = unique(allEdges, 'rows');

        % cell-to-edge connectivity
        c2e = reshape(ic, nCells, nEdgesLoc);
      else
        error("EFEM implemented only for hexa for now");
      end


    end




    function H = computeHeaviside(obj,i,coord)

      f = obj.fractureMesh.surfaces;
      dist = f.center(i,:) - coord;
      n = f.normal(i,:)';
      dn = dist*n;
      assert(all(abs(dn)>1e-10),"Defined fracture passes exactly trough a node. This is not handled yet")
      H = double(dn > 0);

    end

    function Bw = computeCompatibilityMatrix(obj,f,i,coord,gradN)

      H = computeHeaviside(obj,i,coord);
      v = sum(gradN.*H',2);
      n = f.normal(i,:);
      m1 = f.tang1(i,:);
      m2 = f.tang2(i,:);

      v = permute(v,[2 1 3]);

      sym_n_dyad_v = obj.sym_AiBj_plus_AjBi(n,v);
      sym_m1_dyad_v = obj.sym_AiBj_plus_AjBi(m1,v);
      sym_m2_dyad_v = obj.sym_AiBj_plus_AjBi(m2,v);

      Bw = [sym_n_dyad_v, sym_m1_dyad_v, sym_m2_dyad_v];

    end



    function E = computeEquilibriumOperator(obj,f,i)

      n = f.normal(i,:);
      m1 = f.tang1(i,:)';
      m2 = f.tang2(i,:)';

      sym_n_dyad_n = obj.sym_AiBj_plus_AjBi(n,n);
      sym_m1_dyad_n = obj.sym_AiBj_plus_AjBi(m1,n);
      sym_m2_dyad_n = obj.sym_AiBj_plus_AjBi(m2,n);

      A = f.area(i);
      V = obj.grid.cells.volume(f.cutCells(i));

      E = (A/V) * [sym_n_dyad_n, sym_m1_dyad_n, sym_m2_dyad_n];

    end



    function x = localSolve(obj,A,b) %#ok<INUSL>
      % Robust local solve for small condensed systems.
      % Uses a direct solve when possible and falls back to a pseudo-inverse
      % for nearly singular tangential slip tangents.
      if isempty(A)
        x = zeros(0,size(b,2));
        return
      end
      if rcond(A) > 1e-12
        x = A\b;
      else
        x = pinv(A)*b;
      end
    end


    function dtdg = computeDerTractionGap(obj,i,slip,t) %#ok<INUSL>

      % Derivative helper retained for compatibility.  The assembled slip
      % branch uses the augmented direction derivative returned by
      % updateTraction, because q = t_T^proj + cT*Delta g_T.

      state = obj.activeSet.curr(i);
      dtdg = zeros(3,3);
      fId = obj.fractureMesh.surfaces.fracId(i);
      phiF = obj.phi(fId);
      cF = obj.cohesion(fId);
      tauLim = abs(cF - t(1)*tan(phiF));

      switch state
        case {ContactMode.open,ContactMode.stick}
          return
        case {ContactMode.slip,ContactMode.newSlip}
          if norm(slip) > obj.activeSet.tol.sliding
            slipNorm = norm(slip);
            dtdg([2 3],[2 3]) = tauLim * ...
              (slipNorm^2*eye(2) - slip * slip')/slipNorm^3;
          end
      end

    end

  end

  methods (Static)


    function vSym = sym_AiBj_plus_AjBi(a,b)
      % compute symmetric dyadic product of two 3x3 tensor
      % return result into a 6x1 voigt array

      a = reshape(a,3,1,[]);
      b = reshape(b,3,1,[]);

      v = pagemtimes(a,'none',b,'ctranspose') + pagemtimes(b,'none',a,'ctranspose');
      vSym = reshape(v,9,1,[]);
      vSym = vSym([1 5 9 8 7 4],:,:);
      vSym(1:3,:,:) = 0.5*vSym(1:3,:,:);

    end

    function out = getField()
      out = [Poromechanics.getField(), "fractureJump"];
    end

  end

end

