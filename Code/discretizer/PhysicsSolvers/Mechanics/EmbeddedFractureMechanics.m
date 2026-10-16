classdef EmbeddedFractureMechanics < PhysicsSolver

  % Unified two-field EFEM(0) contact solver.
  % isSmooth=true: fixed contact branch during Newton, outer active-set loop.
  % isSmooth=false: branch selected inside each Newton assembly.

  properties
    activeSet = struct("curr",[],"prev",[])
    phi
    cohesion
    fractureMesh
    areaTol = 1e-6
    bcTraction
    mechSolver
    isSmooth = true
    contactAugmentation
  end

  properties (Access = private)
    fldMech
    fldFrac
    trialTraction % augmented trial q, cached by the latest assembly
  end

  methods (Access = public)

    function obj = EmbeddedFractureMechanics(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,varargin)
      obj.mechSolver = Poromechanics(obj.domain);
      obj.mechSolver.registerSolver(varargin{:});
      default = struct('Fracture',struct.empty, ...
        'ActiveSet',missing, ...
        'isSmooth',true, ...
        'augmentationParameter',1.0, ...
        'augmentationNormal',missing, ...
        'augmentationTangential',missing);
      params = readInput(default,varargin{:});
      obj.isSmooth = logical(params.isSmooth);

      if ~ismissing(params.augmentationNormal)
        obj.contactAugmentation(1) = params.augmentationNormal;
      else
        obj.contactAugmentation(1) = params.augmentationParameter;
      end

      if ~ismissing(params.augmentationTangential)
        obj.contactAugmentation(2) = params.augmentationTangential;
      else
        obj.contactAugmentation(2) = params.augmentationParameter;
      end

      defineFractures(obj,params.Fracture);
      nCutCells = obj.fractureMesh.surfaces.num;
      dofm = obj.domain.dofm;
      dofm.registerVariable("fractureJump",entityField.cell,3,"nEntities",nCutCells);
      flds = obj.getField;
      obj.fldMech = dofm.getVariableId(flds(1));
      obj.fldFrac = dofm.getVariableId(flds(2));
      initState(obj);
      initializeActiveSet(obj,nCutCells,params.ActiveSet);
    end

    function assembleSystem(obj,dt)
      assembleSystemEFEM(obj,dt);
    end

    function initialize(obj)
      obj.mechSolver.initialize();
      tIni = computeInitialTraction(obj);
      t = getState(obj,"traction") + tIni(:);
      setState(obj,t,"traction");
      setStateOld(obj,t,"traction");
      setStateInit(obj,t,"traction");
      obj.trialTraction = t;
      obj.activeSet.curr(obj.hasFractureTractionBC()) = ContactMode.open;
      obj.activeSet.prev = obj.activeSet.curr;
    end

    function trac = computeInitialTraction(obj)
      avgStress = getState(obj.mechSolver,"avgStress");
      frac = obj.fractureMesh.surfaces;
      idx = [1;6;5;6;2;4;5;4;3];
      sigma = zeros(3);
      trac = zeros(3*frac.num,1);

      % project average stress onto the embedded fracture
      for i = 1:frac.num
        R = reshape(frac.rotationMatrices(i,:),3,3);
        cellId = frac.cutCells(i);
        sigma(:) = avgStress(cellId,idx);
        t = sigma*frac.normal(i,:)';
        trac(DoFManager.dofExpand(i,3)) = R'*t;
      end
    end

    function timeStepSetup(obj)

      if obj.isSmooth
        hasBC = obj.hasFractureTractionBC();
        closed = obj.activeSet.curr ~= ContactMode.open & ~hasBC;
        obj.activeSet.curr(closed) = ContactMode.stick;

        if isfield(obj.activeSet,'stateChange')
          obj.activeSet.stateChange(:) = 0;
        end
      end

      % set to open fracture element that are initially pressurized
      obj.activeSet.curr(obj.hasFractureTractionBC()) = ContactMode.open;
    end

    function assembleSystemEFEM(obj,dt)

      % allocate
      dofm = obj.domain.dofm;
      frac = obj.fractureMesh.surfaces;
      cells = obj.grid.cells;
      c2f = obj.fractureMesh.cells.cell2fracId;
      subCells = dofm.getFieldCells(obj.fldMech);

      if any(~ismember(frac.cutCells,subCells))
        error('EFEMContact:inactiveCutCell',...
          'All fracture cut cells must belong to the mechanical solver region.');
      end

      nDim = obj.grid.nDim;
      n = sum(nDim^2 * (cells.numVerts(subCells)).^2);
      n1 = nDim^2 * sum(cells.numVerts(frac.cutCells));
      nFracPerCell = cellfun(@numel,c2f);
      n2 = nDim^2 * sum(nFracPerCell.^2,"all");
      nDofU = dofm.getNumbDoF(Poromechanics.getField());
      nDofW = dofm.getNumbDoF(obj.fldFrac);
      asbKuu = assembler(n,nDofU,nDofU);
      asbKuw = assembler(n1,nDofU,nDofW);
      asbKwu = assembler(n1,nDofW,nDofU);
      asbKww = assembler(n2,nDofW,nDofW);
      rhsU = zeros(nDofU,1);
      rhsW = zeros(nDofW,1);

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
      gpMap = obj.domain.gpMap;
      Caug = diag([ ...
        obj.contactAugmentation(1),...
        obj.contactAugmentation(2),...
        obj.contactAugmentation(2)]);

      for cTag = 1:cells.nTag

        % extract constitutive law
        constLaw = obj.domain.materials.getConstitutiveLaw(cTag);

        % extract cells belonging to subregion
        subRegionCells = subCells(cells.tag(subCells) == cTag);

        for vtkId = cells.vtkTypes

          % extract cells of the subregion of homogeneous vtk type
          cellId = cells.VTKType(subRegionCells) == vtkId;
          subCellsLoc = subRegionCells(cellId);

          if isempty(subCellsLoc)
            continue
          end

          cellList = find(cellId);
          elem = FiniteElementType.create( ...
            vtkId,obj.grid,obj.mechSolver.getGaussOrder);

          % get node topology for given vtk type
          topol = obj.grid.getCellNodes(subCellsLoc);
          nG = elem.getNumbGaussPts;

          for i = 1:numel(subCellsLoc)
            el = subCellsLoc(i);
            l = gpMap(el,1);
            fracSurfs = c2f{el};
            nF = numel(fracSurfs);

            % if the cells is not cut by any fracture, nF == 0
            nodes = topol(i,:);
            uDof = dofm.getLocalDoF(obj.fldMech,nodes);
            coords = coordinates(nodes,:);

            % compute strain from continuous displacement field
            [gradN,dJw] = getDerBasisFAndDet(elem,coords);
            B = elem.getStrainMatrix(gradN);
            s.strain(l:l+nG-1,:) = ...
              reshape(pagemtimes(B,du(uDof)),6,nG)';

            % store fracture operators
            Bw = zeros(6,3,nG,nF);
            E = zeros(6,3,1,nF);
            dofW = zeros(3,nF);
            enhancedStrain = zeros(nG,6);

            for fA = 1:nF % loop over all fractures intersecting the cell
              fId = fracSurfs(fA);
              dofA = dofm.getLocalDoF(obj.fldFrac,fId);
              dofW(:,fA) = dofA;

              % precompute enhanced strain matrix and projection operators
              BwLoc = computeCompatibilityMatrix( ...
                obj,frac,fId,coords,gradN);
              Bw(:,:,:,fA) = BwLoc;
              E(:,:,1,fA) = computeEquilibriumOperator(obj,frac,fId);
              enhancedStrain = enhancedStrain + ...
                reshape(pagemtimes(BwLoc,dj(dofA)),6,nG)';
            end

            % total enhanced strain
            s.strain(l:l+nG-1,:) = ...
              s.strain(l:l+nG-1,:) + enhancedStrain;

            % constitutive update (uses the enhanced strain!)
            [sigma,D] = constLaw.constitutiveUpdate( ...
              cellList(i),...
              sOld.stress(l:l+nG-1,:),...
              s.strain(l:l+nG-1,:),...
              dt,...
              time);

            % update stress map
            s.stress(l:l+nG-1,:) = sigma;

            % assemble internal forces for u
            dsigma = sigma - iniStress(l:l+nG-1,:);
            dsigma = reshape(dsigma',6,1,nG);
            fTmp = pagemtimes(B,'ctranspose',dsigma,'none');
            fTmp = fTmp.*reshape(dJw,1,1,[]);
            rhsU(uDof) = rhsU(uDof) + sum(fTmp,3);
            KLoc = obj.mechSolver.computeKloc(B,D,B,dJw);
            asbKuu.localAssembly(uDof,uDof,KLoc);

            % fracture contributions
            for fA = 1:nF
              fId = fracSurfs(fA);
              dofA = dofW(:,fA);
              BwA = Bw(:,:,:,fA);
              EA = E(:,:,:,fA);
              KuwLoc = Poromechanics.computeKloc(B,D,BwA,dJw);
              asbKuw.localAssembly(uDof,dofA,KuwLoc);
              fTmp = pagemtimes(EA,'ctranspose',dsigma,'none');
              fTmp = fTmp.*reshape(dJw,1,1,[]);
              rSigma = sum(fTmp,3);
              A = frac.area(fId);
              projectedTrac = ...
                iniTraction(dofA) + ...
                rSigma/A - ...
                obj.bcTraction(dofA);
              gAug = [jump(dofA(1)); dj(dofA(2:3))];
              q = projectedTrac + Caug*gAug;
              obj.trialTraction(dofA) = q;
              [tracNew,Jproj,mode] = evaluateContactLaw(obj,q,fId);
              s.traction(dofA) = tracNew;
              obj.activeSet.curr(fId) = mode;
              rT = ...
                (tracNew - iniTraction(dofA))*A;

              % prescribed traction on the fracture (pressurization)
              rBC = ...
                obj.bcTraction(dofA)*A;
              rhsW(dofA) = rhsW(dofA) + ...
                rSigma - rT - rBC;
              IminusJ = eye(3) - Jproj;
              KwuLoc = ...
                Poromechanics.computeKloc(EA,D,B,dJw);
              asbKwu.localAssembly( ...
                dofA,uDof,IminusJ*KwuLoc);

              % nested loop (account for fracture intersections)
              for fB = 1:nF
                BwB = Bw(:,:,:,fB);
                dofB = dofW(:,fB);
                KwwLoc = ...
                  obj.mechSolver.computeKloc(EA,D,BwB,dJw);
                KAB = IminusJ*KwwLoc;

                % Augmentation acts only on fracture A's own jump.
                if fA == fB
                  KAB = KAB - A*Jproj*Caug;
                end

                asbKww.localAssembly(dofA,dofB,KAB);
              end % end loop over other fractures in the cell
            end % end loop over fractures in the cell
          end % end loop over cells
        end % end loop over cell shapes
      end % end loop over cell tag

      % update modified state
      setState(obj,s);

      % populate rhs and Jacobian
      obj.domain.J{obj.fldMech,obj.fldMech} = ...
        asbKuu.sparseAssembly;
      obj.domain.J{obj.fldMech,obj.fldFrac} = ...
        asbKuw.sparseAssembly;
      obj.domain.J{obj.fldFrac,obj.fldMech} = ...
        asbKwu.sparseAssembly;
      obj.domain.J{obj.fldFrac,obj.fldFrac} = ...
        asbKww.sparseAssembly;
      obj.domain.rhs{obj.fldMech} = rhsU;
      obj.domain.rhs{obj.fldFrac} = rhsW;
    end

    function initState(obj)
      nCutCells = obj.fractureMesh.surfaces.num;
      state = getState(obj);
      state.fractureJump = zeros(3*nCutCells,1);
      state.traction = zeros(3*nCutCells,1);
      setState(obj,state);
      obj.bcTraction = zeros(3*nCutCells,1);
      obj.activeSet.curr = repmat(ContactMode.stick,nCutCells,1);
      obj.activeSet.prev = obj.activeSet.curr;
      obj.trialTraction = zeros(3*nCutCells,1);
    end

    function hasConfigurationChanged = updateConfiguration(obj)

      if ~obj.isSmooth
        hasConfigurationChanged = false;
        return
      end

      hasConfigurationChanged = updateOuterActiveSet(obj);
    end

    function isReset = resetConfiguration(obj)

      if ~obj.isSmooth
        isReset = false;
        return
      end

      hasBC = obj.hasFractureTractionBC();
      closed = obj.activeSet.curr ~= ContactMode.open & ~hasBC;
      obj.activeSet.curr(closed) = ContactMode.stick;
      obj.activeSet.curr(hasBC) = ContactMode.open;
      isReset = true;
    end

    function advanceState(obj)
      obj.mechSolver.advanceState();
      obj.activeSet.prev = obj.activeSet.curr;
    end

    function order = getGaussOrder(obj)
      order = obj.mechSolver.getGaussOrder;
    end

    function goBackState(obj)
      obj.mechSolver.goBackState();
      obj.activeSet.curr = obj.activeSet.prev;

      if obj.isSmooth && isfield(obj.activeSet,'resetActiveSet') && ...
          obj.activeSet.resetActiveSet
        resetConfiguration(obj);
      end

      obj.activeSet.curr(obj.hasFractureTractionBC()) = ContactMode.open;

      % Next assembly recomputes q from restored stress and jump.
      obj.trialTraction = getStateOld(obj,'traction');
    end

    function updateState(obj,solution)
      obj.mechSolver.updateState(solution);

      if nargin > 1
        dofm = obj.domain.dofm;
        ents = dofm.getActiveEntities(obj.fldFrac,1);
        stateCurr = obj.getState();
        dw = solution(getDoF(dofm,obj.fldFrac));
        stateCurr.fractureJump(ents) = stateCurr.fractureJump(ents) + dw;
        setState(obj,stateCurr);
        end

      if gresLog().getVerbosity >= 2
        nStick = sum(obj.activeSet.curr == ContactMode.stick);
        nSlip = sum(obj.activeSet.curr == ContactMode.slip | ...
          obj.activeSet.curr == ContactMode.newSlip);
        nOpen = sum(obj.activeSet.curr == ContactMode.open);
        fprintf('%s: active set ',class(obj));
        fprintf('stick = %i, slip = %i, open = %i\n',nStick,nSlip,nOpen);
      end
    end

    function applyBC(obj,bcId,t)
      obj.mechSolver.applyBC(bcId,t);
    end

    function applyDirVal(obj,bcId,t)
      obj.mechSolver.applyDirVal(bcId,t);
    end

    function out = isLinear(obj)
      out = obj.isSmooth && obj.mechSolver.isLinear() && ...
        all(obj.activeSet.curr == ContactMode.stick);
    end

    function out = isSymmetric(obj) %#ok<MANU>

      % In stick K_wu=0 while B_1=K_uw is retained: nonsymmetric.
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
      cellDataFrac = repmat(struct('name',1,'data',1),1,1);
      cellDataFrac(1).name = 'isCellFractured';
      isCutCell = false(obj.grid.cells.num,1);
      isCutCell(obj.fractureMesh.surfaces.cutCells) = true;
      cellDataFrac(1).data = double(isCutCell);
      cellData = OutState.mergeOutFields(cellData,cellDataFrac);
      s = obj.domain.state.interpolate(fac);
      jump = s.fractureJump;
      jumpDiff = s.fractureJump-obj.domain.getStateOld.fractureJump;
      trac = s.traction(:);
      cellStr = repmat(struct('name',1,'data',1),4,1);
      cellStr(1).name='fractureJump';
      cellStr(1).data=[jump(1:3:end),jump(2:3:end),jump(3:3:end)];
      cellStr(2).name='traction';
      cellStr(2).data=[trac(1:3:end),trac(2:3:end),trac(3:3:end)];
      cellStr(3).name='fractureState';
      cellStr(3).data=double(obj.activeSet.curr);
      cellStr(4).name='jumpSlip';
      cellStr(4).data=[jumpDiff(2:3:end),jumpDiff(3:3:end)];
      obj.domain.outstate.writeVTKfile(obj.domain.vtmBlock,'EmbeddedFractures', ...
        obj.fractureMesh,time,[],[],[],cellStr);
    end
  end

  methods (Access=private)

    function [traction,J,mode] = evaluateContactLaw(obj,q,fId)

      % Self-contained strategy switch. The assembly always uses the same
      % projected traction and its derivative with respect to q.
      % classify contact according to algorithm strategy (smooth vs
      % semi-smooth)
      % return updated traction
      % derivative of traction with respect to the augmented trial q
      % contact mode
      if obj.hasFractureTractionBC(fId)
        mode = ContactMode.open;
      elseif obj.isSmooth
        mode = obj.activeSet.curr(fId);
      else
        mode = classifyContact(obj,q,obj.fractureMesh.surfaces.fracId(fId));
      end

      [traction,J] = evaluateContactBranch(obj,q,...
        obj.fractureMesh.surfaces.fracId(fId),mode,fId);
    end

    function mode = classifyContact(obj,q,fracId)

      % Same tolerance policy as SolidMechanicsContact's in-Newton
      % classification: strict normal test and relative tangential allowance.
      if q(1) > 0
        mode = ContactMode.open;
        return
      end

      tauLim = max(obj.cohesion(fracId)-tan(obj.phi(fracId))*q(1),0);
      tau = norm(q(2:3));

      % if tau >= tauLim
      %   tau = tau*(1-obj.activeSet.tol.tangentialViolation);
      % end
      if tau <= tauLim
        mode = ContactMode.stick;
      else
        mode = ContactMode.slip;
      end
    end

    function [traction,J] = evaluateContactBranch(obj,q,fracId,mode,fId)

      % compute updated traction and jacobian block dtdg
      traction = zeros(3,1);
      J = zeros(3);

      if mode == ContactMode.open
        return
      elseif mode == ContactMode.stick

        % Fixed-stick trial must retain inadmissible reactions so the outer
        % update can detect slip/open. No clipping or jump overwriting.
        traction = q;
        J = eye(3);
        return
      elseif mode ~= ContactMode.slip && mode ~= ContactMode.newSlip
        error('EFEMContact:mode','Unsupported contact branch.');
      end

      traction(1) = q(1);
      J(1,1) = 1;
      mu = tan(obj.phi(fracId));
      tauRaw = obj.cohesion(fracId)-mu*q(1);

      if obj.isSmooth

        % Fixed-slip continuation matches the old outer-loop formulation.
        % Select the positive-side derivative of abs at tauRaw = 0.
        tauLim = abs(tauRaw);
        dTauDqN = -mu;

        if tauRaw < 0
          dTauDqN = mu;
        end

      else

        % Semi-smooth opening is selected by classifyContact; retain the
        % clamped Coulomb bound and its corresponding derivative.
        tauLim = max(tauRaw,0);
        dTauDqN = -mu*double(tauRaw > 0);
      end

      qT = q(2:3);
      qNorm = norm(qT);

      if qNorm > obj.activeSet.tol.sliding
        direction = qT/qNorm;
        traction(2:3) = tauLim*direction;
        J(2:3,1) = dTauDqN*direction;
        J(2:3,2:3) = (tauLim/qNorm)*(eye(2)-direction*direction');

      elseif obj.isSmooth

        % A fixed slip branch can encounter a vanishing augmented direction
        % during unloading/reversal. Retain a reaction direction, as in the
        % old EFEMaugmented solver, instead of dropping friction to zero.
        direction = getFixedSlipDirection(obj,qT,fId);
        traction(2:3) = tauLim*direction;
        J(2:3,1) = dTauDqN*direction;

        % The fallback direction is held fixed in this local linearization.
        J(2:3,2:3) = zeros(2);
      end

      % The semismooth cutoff behavior is intentionally preserved.
    end

    function direction = getFixedSlipDirection(obj,qT,fId)
      id = DoFManager.dofExpand(fId,3);
      state = getState(obj);
      stateOld = getStateOld(obj);
      tCurr = state.traction(id(2:3));
      tOld = stateOld.traction(id(2:3));

      if norm(tCurr) > 0
        direction = tCurr/norm(tCurr);
      elseif norm(tOld) > 0
        direction = tOld/norm(tOld);
      elseif norm(qT) > 0
        direction = qT/norm(qT);
      else
        direction = [1;0];
      end
    end

    function changed = updateOuterActiveSet(obj)
      oldModes = obj.activeSet.curr;
      newModes = oldModes;
      frac = obj.fractureMesh.surfaces;
      hasBC = obj.hasFractureTractionBC();
      jump = obj.getState.fractureJump;

      % The driver assembles the final iterate before checking convergence.
      % Use its unprojected augmented trial: a returned slip traction lies
      % on the friction bound and hides unloading back into stick.
      trial = obj.trialTraction;

      for i = 1:frac.num

        if hasBC(i)
          newModes(i) = ContactMode.open;
        else
          id = DoFManager.dofExpand(i,3);
          f = frac.fracId(i);
          q = trial(id);
          tauLim = abs(obj.cohesion(f)-tan(obj.phi(f))*q(1));
          gN = jump(id(1));
          newModes(i) = updateContactState(oldModes(i),q,tauLim,gN,obj.activeSet.tol);

          if newModes(i) == ContactMode.slip && oldModes(i) == ContactMode.stick
            newModes(i) = ContactMode.newSlip;
          end
        end
      end

      oldSliding = oldModes == ContactMode.slip | oldModes == ContactMode.newSlip;
      newSliding = newModes == ContactMode.slip | newModes == ContactMode.newSlip;
      changedFaces = oldModes ~= newModes & ~(oldSliding & newSliding);
      obj.activeSet.curr = newModes;

      % Counters are diagnostic only. Do not suppress an inadmissible branch
      % change because an area/count tolerance has been exceeded.
      if isfield(obj.activeSet,'stateChange')
        obj.activeSet.stateChange(changedFaces) = ...
          obj.activeSet.stateChange(changedFaces)+1;
      end

      changed = any(changedFaces);
      gresLog().log(3,'%s: outer active set changed on %i faces.\n',...
        class(obj),sum(changedFaces));
    end

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
      f.fracId        = f.fracId(id);
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
      c = fMesh.cells;
      nC = obj.grid.cells.num;
      cell2fracId = cell(nC,1);

      for i = 1:nC
        cell2fracId{i} = reshape(find(f.cutCells == i),1,[]);
      end

      c.cell2fracId = cell2fracId;
      fMesh.cells = c;
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

    % HELPER FUNCTIONS

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
