classdef EFEMsemismooth < PhysicsSolver

  % Monolithic semismooth EFEM(0) contact solver.
  %
  % All fracture-jump components are solved in one Newton system. Contact
  % modes are obtained implicitly from a semismooth projection and are only
  % stored for diagnostics/output; no outer active-set loop is performed.

  properties
    activeSet = struct("curr",[],"prev",[]) % diagnostic only
    phi
    cohesion
    fractureMesh
    areaTol = 1e-6
    bcTraction
    mechSolver
    fixActiveSet = true
    augmentation_n = 1e3
    augmentation_t = 1e3
    projectionTol = 1e-12
  end

  properties (Access = private)
    fldMech
    fldFrac
  end

  methods (Access = public)

    function obj = EFEMsemismooth(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,varargin)
      obj.mechSolver = Poromechanics(obj.domain);
      obj.mechSolver.registerSolver(varargin{:});

      default = struct('Fracture',struct.empty, ...
        'augmentationNormal',1e3, ...
        'augmentationTangential',1e3, ...
        'projectionTolerance',1e-12);
      params = readInput(default,varargin{:});
      obj.augmentation_n = params.augmentationNormal;
      obj.augmentation_t = params.augmentationTangential;
      obj.projectionTol = params.projectionTolerance;

      defineFractures(obj,params.Fracture);
      nCutCells = obj.fractureMesh.surfaces.num;
      dofm = obj.domain.dofm;
      dofm.registerVariable("fractureJump",entityField.cell,3,"nEntities",nCutCells);

      flds = obj.getField;
      obj.fldMech = dofm.getVariableId(flds(1));
      obj.fldFrac = dofm.getVariableId(flds(2));
      initState(obj);
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
    end

    function trac = computeInitialTraction(obj)
      avgStress = getState(obj.mechSolver,"avgStress");
      frac = obj.fractureMesh.surfaces;
      idx = [1;6;5;6;2;4;5;4;3];
      sigma = zeros(3);
      trac = zeros(3*frac.num,1);
      for i = 1:frac.num
        R = reshape(frac.rotationMatrices(i,:),3,3);
        cellId = frac.cutCells(i);
        sigma(:) = avgStress(cellId,idx);
        t = sigma*frac.normal(i,:)';
        trac(DoFManager.dofExpand(i,3)) = R'*t;
      end
    end

    function timeStepSetup(obj) %#ok<MANU>
      % No contact configuration is reset: the semismooth map is evaluated
      % directly at every Newton iterate.
    end

    function assembleSystemEFEM(obj,dt)
      % Assemble semi-smooth EFEM system.
      %
      % Supports:
      %   - multiple bulk element types;
      %   - multiple material regions;
      %   - multiple fracture surfaces intersecting the same bulk cell;
      %   - full fracture-fracture coupling through the bulk stress.

      % allocate
      dofm = obj.domain.dofm;
      frac = obj.fractureMesh.surfaces;
      cells = obj.grid.cells;

      c2f = obj.fractureMesh.cells.cell2fracId;

      subCells = dofm.getFieldCells(obj.fldMech);

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
        obj.augmentation_n,...
        obj.augmentation_t,...
        obj.augmentation_t]);

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

            for fA = 1:nF

              fId = fracSurfs(fA);

              dofA = dofm.getLocalDoF(obj.fldFrac,fId);
              dofW(:,fA) = dofA;

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

            % constitutive update
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

              % -------------------------------------------------------------
              % Kuw contribution associated with fracture A
              % -------------------------------------------------------------

              KuwLoc = Poromechanics.computeKloc(B,D,BwA,dJw);
              asbKuw.localAssembly(uDof,dofA,KuwLoc);

              % -------------------------------------------------------------
              % projected bulk stress on fracture A
              % -------------------------------------------------------------

              fTmp = pagemtimes(EA,'ctranspose',dsigma,'none');
              fTmp = fTmp.*reshape(dJw,1,1,[]);
              rSigma = sum(fTmp,3);

              A = frac.area(fId);

              projectedTrac = ...
                iniTraction(dofA) + ...
                rSigma/A - ...
                obj.bcTraction(dofA);

              % -------------------------------------------------------------
              % semi-smooth traction projection
              % -------------------------------------------------------------

              if obj.hasFractureTractionBC(fId)

                tracNew = zeros(3,1);
                Jproj = zeros(3);
                mode = ContactMode.open;

              else

                gAug = [ ...
                  jump(dofA(1)); ...
                  dj(dofA(2:3))];

                q = projectedTrac + Caug*gAug;

                [tracNew,Jproj,mode] = ...
                  projectContactTraction( ...
                  obj,q,frac.fracId(fId));

              end

              s.traction(dofA) = tracNew;
              obj.activeSet.curr(fId) = mode;

              % -------------------------------------------------------------
              % fracture residual
              % -------------------------------------------------------------

              rT = ...
                (tracNew - iniTraction(dofA))*A;

              rBC = ...
                obj.bcTraction(dofA)*A;

              rhsW(dofA) = rhsW(dofA) + ...
                rSigma - rT - rBC;

              % -------------------------------------------------------------
              % Kwu
              % -------------------------------------------------------------

              IminusJ = eye(3) - Jproj;

              KwuLoc = ...
                Poromechanics.computeKloc(EA,D,B,dJw);

              asbKwu.localAssembly( ...
                dofA,uDof,IminusJ*KwuLoc);

              % -------------------------------------------------------------
              % Kww
              %
              % Fracture equation A depends on the jumps of every fracture B
              % intersecting the same bulk cell.
              % -------------------------------------------------------------

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

              end

            end

          end

        end

      end

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



    %     function assembleSystemEFEM(obj,dt)
    %   % Assemble EFEM(0) system with state-dependent active fracture DOFs.
    %   % Multiple fracture surfaces can intersect the same bulk cell.
    %
    %   % allocate
    %   dofm = obj.domain.dofm;
    %   frac = obj.fractureMesh.surfaces;
    %   cells = obj.grid.cells;
    %
    %   c2f = obj.fractureMesh.cells.cell2fracId;
    %
    %   subCells = dofm.getFieldCells(obj.fldMech);
    %
    %   nDim = obj.grid.nDim;
    %   n = sum(nDim^2 * (obj.grid.cells.numVerts(subCells)).^2);
    %   n1 = nDim^2 * sum(cells.numVerts(frac.cutCells));
    %   nFracPerCell = cellfun(@numel,c2f);
    %   n2 = nDim^2 * sum(nFracPerCell.^2,"all");
    %   nDofU = dofm.getNumbDoF(Poromechanics.getField());
    %   nDofW = dofm.getNumbDoF(obj.fldFrac);
    %
    %   asbKuu = assembler(n,nDofU,nDofU);
    %   asbKuw = assembler(n1,nDofU,nDofW);
    %   asbKwu = assembler(n1,nDofW,nDofU);
    %   asbKww = assembler(n2,nDofW,nDofW);
    %   rhsU = zeros(nDofU,1);
    %   rhsW = zeros(nDofW,1);
    %
    %   % get state variables
    %   s = getState(obj);
    %   sOld = getStateOld(obj);
    %   time = s.time;
    %   iniStress = getStateInit(obj,'stress');
    %   iniTraction = getStateInit(obj,'traction');
    %   jump = s.fractureJump;
    %   jumpOld = sOld.fractureJump;
    %   du = s.displacements - sOld.displacements;
    %   dj = jump - jumpOld;
    %
    %   coordinates = obj.grid.coordinates;
    %   gpMap = obj.domain.gpMap;
    %
    %   for cTag = 1:cells.nTag
    %
    %     % extract the constitutive law
    %     constLaw = obj.domain.materials.getConstitutiveLaw(cTag);
    %
    %     % extract cells belonging to subregion
    %     subRegionCells = subCells(cells.tag(subCells) == cTag);
    %
    %     for vtkId = cells.vtkTypes
    %
    %       % extract cells of the subregion of homogeneous vtk type
    %       cellId = cells.VTKType(subRegionCells) == vtkId;
    %       subCellsLoc = subRegionCells(cellId);
    %
    %       if isempty(subCellsLoc)
    %         continue
    %       end
    %
    %       cellList = find(cellId);
    %
    %       elem = FiniteElementType.create(vtkId,obj.grid,...
    %         obj.mechSolver.getGaussOrder);
    %
    %       % get node topology for given vtk type
    %       topol = obj.grid.getCellNodes(subCellsLoc);
    %
    %       nG = elem.getNumbGaussPts;
    %
    %       for i = 1:numel(subCellsLoc)
    %
    %         el = subCellsLoc(i);
    %         l = gpMap(el,1);
    %
    %         fracSurfs = c2f{el};
    %         nF = numel(fracSurfs);
    %
    %         nodes = topol(i,:);
    %         uDof = dofm.getLocalDoF(obj.fldMech,nodes);
    %         coords = coordinates(nodes,:);
    %
    %         % compute strain from continuous displacement field
    %         [gradN,dJw] = getDerBasisFAndDet(elem,coords);
    %
    %         B = elem.getStrainMatrix(gradN);
    %         s.strain(l:l+nG-1,:) = ...
    %           reshape(pagemtimes(B,du(uDof)),6,nG)';
    %
    %         enhancedStrain = zeros(nG,6);
    %
    %         Bw = zeros(6,3,nG,nF);
    %         E = zeros(6,3,1,nF);
    %         dofW = zeros(3,nF);
    %         fracState = zeros(nF,1);
    %
    %         for fA = 1:nF
    %
    %           fId = fracSurfs(fA);
    %
    %           % grab fracture DOFs and operators
    %           dofW(:,fA) = dofm.getLocalDoF(obj.fldFrac,fId);
    %           BwLoc = computeCompatibilityMatrix(obj,frac,fId,coords,gradN);
    %           Bw(:,:,:,fA) = BwLoc;
    %           E(:,:,1,fA) = computeEquilibriumOperator(obj,frac,fId);
    %
    %           state = obj.activeSet.curr(fId);
    %           if obj.hasFractureTractionBC(fId)
    %             state = ContactMode.open;
    %             obj.activeSet.curr(fId) = ContactMode.open;
    %           end
    %           fracState(fA) = state;
    %
    %           % Enhance strain only with jump components allowed by the
    %           % current contact state.
    %           dofA = dofW(:,fA);
    %           djEff = dj(dofA);
    %
    %           if state == ContactMode.stick
    %             djEff(:) = 0;
    %             s.fractureJump(dofA(1)) = 0;
    %             s.fractureJump(dofA(2:3)) = jumpOld(dofA(2:3));
    %           end
    %
    %           enhancedStrain = enhancedStrain + ...
    %             reshape(pagemtimes(BwLoc,djEff),6,nG)';
    %
    %         end
    %
    %         % total enhanced strain
    %         s.strain(l:l+nG-1,:) = ...
    %           s.strain(l:l+nG-1,:) + enhancedStrain;
    %
    %         % constitutive update
    %         [sigma,D] = constLaw.constitutiveUpdate(cellList(i),...
    %           sOld.stress(l:l+nG-1,:),...
    %           s.strain(l:l+nG-1,:),...
    %           dt,...
    %           time);
    %
    %         % update stress map
    %         s.stress(l:l+nG-1,:) = sigma;
    %
    %         % assemble internal forces for u
    %         dsigma = sigma - iniStress(l:l+nG-1,:);
    %         dsigma = reshape(dsigma',6,1,nG);
    %         fTmp = pagemtimes(B,'ctranspose',dsigma,'none');
    %         fTmp = fTmp.*reshape(dJw,1,1,[]);
    %         fLoc = sum(fTmp,3);
    %         rhsU(uDof) = rhsU(uDof) + fLoc;
    %
    %         KLoc = obj.mechSolver.computeKloc(B,D,B,dJw);
    %         asbKuu.localAssembly(uDof,uDof,KLoc);
    %
    %         % Mechanical coupling columns depend on the state of each
    %         % fracture contributing to the enhanced strain.
    %         for fB = 1:nF
    %
    %           stateB = fracState(fB);
    %           if stateB == ContactMode.stick
    %             continue
    %           end
    %
    %           BwB = Bw(:,:,:,fB);
    %           dofB = dofW(:,fB);
    %           KuwLoc = Poromechanics.computeKloc(B,D,BwB,dJw);
    %
    %           if stateB == ContactMode.slip || ...
    %               stateB == ContactMode.newSlip
    %             asbKuw.localAssembly(uDof,dofB(2:3),KuwLoc(:,2:3));
    %           elseif stateB == ContactMode.open
    %             asbKuw.localAssembly(uDof,dofB,KuwLoc);
    %           else
    %             error('Unsupported contact state in EFEM assembly')
    %           end
    %
    %         end
    %
    %         for fA = 1:nF
    %
    %           fId = fracSurfs(fA);
    %           dofA = dofW(:,fA);
    %           stateA = fracState(fA);
    %           EA = E(:,:,:,fA);
    %
    %           % projected stress residual on fracture A
    %           fTmp = pagemtimes(EA,'ctranspose',dsigma,'none');
    %           fTmp = fTmp.*reshape(dJw,1,1,[]);
    %           rSigma = sum(fTmp,3);
    %
    %           tCurr = s.traction(dofA);
    %           tOld = sOld.traction(dofA);
    %           fracId = frac.fracId(fId);
    %
    %           if stateA == ContactMode.stick
    %
    %             % No fracture jump is solved in stick. The traction is the
    %             % reaction obtained from the projected bulk stress.
    %             tracNew = iniTraction(dofA) + ...
    %               rSigma/frac.area(fId) - obj.bcTraction(dofA);
    %             s.traction(dofA) = tracNew;
    %
    %             asbKww.localAssembly(dofA,dofA,eye(3));
    %
    %             jStick = [jump(dofA(1)); dj(dofA(2:3))];
    %             rhsW(dofA) = rhsW(dofA) + jStick;
    %
    %           elseif stateA == ContactMode.slip || ...
    %               stateA == ContactMode.newSlip
    %
    %             projectedTrac = iniTraction(dofA) + ...
    %               rSigma/frac.area(fId) - obj.bcTraction(dofA);
    %
    %             [tracNew,slipDir,dDir,tauLim,tauSign] = ...
    %               updateTraction(obj,fId,fracId,projectedTrac,tOld,tCurr,...
    %               s.fractureJump(dofA),jumpOld(dofA));
    %             s.traction(dofA) = tracNew;
    %
    %             rT = (tracNew - iniTraction(dofA))*frac.area(fId);
    %             rBC = obj.bcTraction(dofA)*frac.area(fId);
    %             rhsLoc = rSigma - rT - rBC;
    %
    %             idN = 1;
    %             idT = [2 3];
    %             Afrac = frac.area(fId);
    %             tanPhi = tan(obj.phi(fracId));
    %             Iaug = eye(2) - tauLim*dDir;
    %
    %             KwuLoc = Poromechanics.computeKloc(EA,D,B,dJw);
    %             KTu = Iaug*KwuLoc(idT,:) + ...
    %               tanPhi*tauSign*(slipDir*KwuLoc(idN,:));
    %
    %             asbKwu.localAssembly(dofA(idT),uDof,KTu);
    %             rhsW(dofA(idT)) = rhsW(dofA(idT)) + rhsLoc(idT);
    %
    %             % Normal contact closure: w_N = 0.
    %             asbKww.localAssembly(dofA(idN),dofA(idN),1);
    %             rhsW(dofA(idN)) = rhsW(dofA(idN)) + ...
    %               s.fractureJump(dofA(idN));
    %
    %             % All fracture pairs contribute through the common stress.
    %             for fB = 1:nF
    %
    %               stateB = fracState(fB);
    %               if stateB == ContactMode.stick
    %                 continue
    %               end
    %
    %               BwB = Bw(:,:,:,fB);
    %               dofB = dofW(:,fB);
    %               KwwLoc = obj.mechSolver.computeKloc(EA,D,BwB,dJw);
    %
    %               if stateB == ContactMode.slip || ...
    %                   stateB == ContactMode.newSlip
    %                 colsB = idT;
    %               elseif stateB == ContactMode.open
    %                 colsB = 1:3;
    %               else
    %                 error('Unsupported contact state in EFEM assembly')
    %               end
    %
    %               KTB = Iaug*KwwLoc(idT,colsB) + ...
    %                 tanPhi*tauSign*(slipDir*KwwLoc(idN,colsB));
    %
    %               if fA == fB
    %                 idTang = find(colsB == 2 | colsB == 3);
    %                 KTB(:,idTang) = KTB(:,idTang) - ...
    %                   Afrac*tauLim*obj.augmentation_t*dDir(:,colsB(idTang)-1);
    %               end
    %
    %               asbKww.localAssembly(dofA(idT),dofB(colsB),KTB);
    %
    %             end
    %
    %           elseif stateA == ContactMode.open
    %
    %             % Open: all jump components are global and projected
    %             % traction is zero.
    %             tracNew = zeros(3,1);
    %             s.traction(dofA) = tracNew;
    %
    %             KwuLoc = Poromechanics.computeKloc(EA,D,B,dJw);
    %             asbKwu.localAssembly(dofA,uDof,KwuLoc);
    %
    %             rT = (tracNew - iniTraction(dofA))*frac.area(fId);
    %             rBC = obj.bcTraction(dofA)*frac.area(fId);
    %             rhsLoc = rSigma - rT - rBC;
    %             rhsW(dofA) = rhsW(dofA) + rhsLoc;
    %
    %             for fB = 1:nF
    %
    %               stateB = fracState(fB);
    %               if stateB == ContactMode.stick
    %                 continue
    %               end
    %
    %               BwB = Bw(:,:,:,fB);
    %               dofB = dofW(:,fB);
    %               KwwLoc = obj.mechSolver.computeKloc(EA,D,BwB,dJw);
    %
    %               if stateB == ContactMode.slip || ...
    %                   stateB == ContactMode.newSlip
    %                 colsB = [2 3];
    %               elseif stateB == ContactMode.open
    %                 colsB = 1:3;
    %               else
    %                 error('Unsupported contact state in EFEM assembly')
    %               end
    %
    %               asbKww.localAssembly(dofA,dofB(colsB),...
    %                 KwwLoc(:,colsB));
    %
    %             end
    %
    %           else
    %             error('Unsupported contact state in EFEM assembly')
    %           end
    %
    %         end
    %
    %       end
    %
    %     end
    %
    %   end
    %
    %   % update modified state
    %   setState(obj,s);
    %
    %   % populate rhs and jacobian
    %   obj.domain.J{obj.fldMech,obj.fldMech} = asbKuu.sparseAssembly;
    %   obj.domain.J{obj.fldMech,obj.fldFrac} = asbKuw.sparseAssembly;
    %   obj.domain.J{obj.fldFrac,obj.fldMech} = asbKwu.sparseAssembly;
    %   obj.domain.J{obj.fldFrac,obj.fldFrac} = asbKww.sparseAssembly;
    %   obj.domain.rhs{obj.fldMech} = rhsU;
    %   obj.domain.rhs{obj.fldFrac} = rhsW;
    %
    % end

    function initState(obj)
      nCutCells = obj.fractureMesh.surfaces.num;
      state = getState(obj);
      state.fractureJump = zeros(3*nCutCells,1);
      state.traction = zeros(3*nCutCells,1);
      setState(obj,state);
      obj.bcTraction = zeros(3*nCutCells,1);
      obj.activeSet.curr = repmat(ContactMode.stick,nCutCells,1);
      obj.activeSet.prev = obj.activeSet.curr;
    end

    function hasConfigurationChanged = updateConfiguration(obj) %#ok<MANU>
      hasConfigurationChanged = false;
    end

    function isReset = resetConfiguration(obj) %#ok<MANU>
      isReset = false;
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
      out = false;
    end

    function out = isSymmetric(obj)
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

    function [traction,J,mode] = projectContactTraction(obj,q,fracId)
      % Semismooth projection associated with
      %   q_N = t_N^proj + c_N g_N,
      %   q_T = t_T^proj + c_T Delta g_T.
      % Compression is negative. q_N >= 0 identifies an open interface.
      mu = tan(obj.phi(fracId));
      coh = obj.cohesion(fracId);
      J = zeros(3);
      traction = zeros(3,1);

      if q(1) >= obj.projectionTol
        mode = ContactMode.open;
        return
      end

      traction(1) = q(1);
      tauLim = max(coh-mu*q(1),0);
      qT = q(2:3);
      qTnorm = norm(qT);

      if qTnorm <= tauLim + obj.projectionTol
        traction(2:3) = qT;
        J = eye(3);
        mode = ContactMode.stick;
        return
      end

      if qTnorm > obj.projectionTol
        e = qT/qTnorm;
        traction(2:3) = tauLim*e;
        J(1,1) = 1;
        J(2:3,1) = -mu*e;
        J(2:3,2:3) = (tauLim/qTnorm)*(eye(2)-e*e');
      else
        J(1,1) = 1;
      end
      mode = ContactMode.slip;
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

