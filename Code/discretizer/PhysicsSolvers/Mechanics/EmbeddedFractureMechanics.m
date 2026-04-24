classdef EmbeddedFractureMechanics < PhysicsSolver

  % solver for embedded tractions implementing the EFEM(0) formulation
  % Cusini et al (2021).

  properties
    activeSet = struct("curr",[],"prev",[])
    penalty_n               % penalty parameter for normal direction
    penalty_t               % penalty parameter for tangential direction
    phi                     % friction angle in radians for each fracture
    cohesion                % the cohesion of each fracture
    fractureMesh            % a 2D mesh object with cut cell topology
    areaTol = 1e-6;         % minimum area of a fracture element
    bcTraction
    mechSolver              % handle to the Poromechanics solver

  end

  properties (Access = private)
    fldMech
    fldFrac
  end

  methods (Access = public)

    function obj = EmbeddedFractureMechanics(domain)

      % call physicsSolver constructor
      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,varargin)

      obj.mechSolver = Poromechanics(obj.domain);
      obj.mechSolver.registerSolver(varargin{:});


      default = struct('penaltyNormal',1e8,...
        'penaltyTangential',1e8,...
        'Fracture',struct.empty,...
        'ActiveSet',missing);

      params = readInput(default,varargin{:});

      obj.penalty_n = params.penaltyNormal;
      obj.penalty_t = params.penaltyTangential;

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
      setStateInit(obj,t,"traction");

    end

    function trac = computeInitialTraction(obj)
      % initialize traction for cell stress (average)
      avgStress = getState(obj.mechSolver,"avgStress");
      frac = obj.fractureMesh.surfaces;
      idx = [1;6;5;6;2;4;5;4;3];
      sigma = zeros(3);
      trac = zeros(frac.num,1);
      for i = 1:frac.num
        R = reshape(frac.rotationMatrices(i,:),3,3);
        cellId = frac.cutCells(i);
        sigma(:) = avgStress(cellId,idx);
        n = frac.normal(i,:);
        t = sigma*n';
        trac(dofId(i,3)) = R'*t;
      end
    end


    function assembleSystemEFEM(obj,dt)

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

      % get state variables
      t = obj.domain.state.t;
      s = getState(obj);
      sOld = getStateOld(obj);
      iniStress = getStateInit(obj,'stress');
      iniTraction = getStateInit(obj,'traction'); % use this!
      jump = s.fractureJump;
      du = s.displacements - sOld.displacements;
      dj = s.fractureJump - sOld.fractureJump;

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

        % compute strain
        [gradN,dJw] = getDerBasisFAndDet(elem,coords);

        B = elem.getStrainMatrix(gradN);
        s.strain(l:l+nG-1,:) = reshape(pagemtimes(B,du(uDof)),6,nG)';

        if isCellCut

          % grab degree of freedom
          wDof = dofm.getLocalDoF(obj.fldFrac,f);

          % compute Bw matrix (compatibility operator, 6x3)
          Bw = computeCompatibilityMatrix(obj,frac,f,coords,gradN);

          % enhance strain
          enhancedStrain = reshape(pagemtimes(Bw,dj(wDof)),6,nG)';

          % compute E matrix (equilibrium operator, 6x3)
          E = computeEquilibriumOperator(obj,frac,f);

          s.strain(l:l+nG-1,:) = s.strain(l:l+nG-1,:) + enhancedStrain;

        end


        % constitutive update 
        [D, sigma, status] = obj.domain.materials.updateMaterial( ...
          cells.tag(el), ...
          sOld.stress(l:l+nG-1,:), ...
          s.strain(l:l+nG-1,:), ...
          dt, sOld.status(l:l+nG-1,:), el, t);

        % update stress map and gp counter
        s.status(l:l+nG-1,:) = status;
        s.stress(l:(l+nG-1),:) = sigma;

        % assemble internal forces
        dsigma = sigma - iniStress(l:l+nG-1,:);
        dsigma = reshape(dsigma',6,1,nG);
        fTmp = pagemtimes(B,'ctranspose',dsigma,'none');
        fTmp = fTmp.*reshape(dJw,1,1,[]);
        fLoc = sum(fTmp,3);
        rhsU(uDof) = rhsU(uDof)+fLoc;

        KLoc = obj.mechSolver.computeKloc(B,D,B,dJw);
        asbKuu.localAssembly(uDof,uDof,KLoc);

        if isCellCut
          % assemble the efem blocks
          KuwLoc = Poromechanics.computeKloc(B,D,Bw,dJw);
          KwuLoc = Poromechanics.computeKloc(E,D,B,dJw);
          KwwLoc = Poromechanics.computeKloc(E,D,Bw,dJw);

          slip = dj(wDof([2;3])); % why not dj?
          trac = s.traction(wDof);
          dtrac = trac - iniTraction(wDof);
          dtdg = computeDerTractionGap(obj,f,slip,trac(1));

          KwwLoc = KwwLoc - dtdg*frac.area(f);

          % assemble local contributions
          asbKuw.localAssembly(uDof,wDof,KuwLoc);
          asbKwu.localAssembly(wDof,uDof,KwuLoc);
          asbKww.localAssembly(wDof,wDof,KwwLoc);

          % assemble rhsW (use computed stress tensor)
          sigma = reshape(sigma',6,1,nG);
          rT = dtrac*frac.area(f);
          fTmp = pagemtimes(E,'ctranspose',sigma,'none');
          fTmp = fTmp.*reshape(dJw,1,1,[]);
          rSigma = sum(fTmp,3);
          rBC = obj.bcTraction(wDof)*frac.area(f);
          rhsLoc = rSigma - rT - rBC;
          rhsW(wDof) = rhsW(wDof) + rhsLoc;

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
      setState(obj,state);
      obj.bcTraction = zeros(3*nCutCells,1);
      obj.activeSet.curr = repmat(ContactMode.open,nCutCells,1);
      obj.activeSet.prev = obj.activeSet.curr;

    end


    function hasConfigurationChanged = updateConfiguration(obj)

      oldActiveSet = obj.activeSet.curr;

      traction = getState(obj,"traction");
      displacementJump = getState(obj,"fractureJump");

      f = obj.fractureMesh.surfaces;


      for i = 1:numel(obj.activeSet.curr)

        fracId = f.fracId(i);

        state = obj.activeSet.curr(i);

        id = DoFManager.dofExpand(i,3);

        t = traction(id);
        g_n = displacementJump(id(1));

        limitTraction = abs(obj.cohesion(fracId) - tan(obj.phi(fracId))*t(1));

        obj.activeSet.curr(i) = updateContactState(state,t,...
          limitTraction, ...
          g_n,...
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
        fprintf('%s: Active set \n',class(obj));
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

      if hasConfigurationChanged
        updateTraction(obj);
      end

    end



    function isReset = resetConfiguration(obj)

      toReset = obj.activeSet.curr(:) ~= ContactMode.open;
      obj.activeSet.curr(toReset) = ContactMode.stick;

      isReset = true;

    end


    function advanceState(obj)
      % Set converged state to current state after newton convergence
      obj.mechSolver.advanceState()
      obj.activeSet.old = obj.activeSet.curr;
    end



    function goBackState(obj)

      % reset state to beginning of time step
      obj.mechSolver.goBackState();

      obj.activeSet.curr = obj.activeSet.prev;
      obj.NLIter = 0;

      if obj.activeSet.resetActiveSet
        resetConfiguration(obj);
      end

    end

    function updateState(obj,solution)

      % Update state structure with last solution increment
      obj.mechSolver.updateState(solution);

      dofm = obj.domain.dofm;
      ents = dofm.getActiveEntities(obj.fldFrac,1);
      stateCurr = obj.getState();
      %stateOld = obj.getStateOld();

      if nargin > 1
        % apply newton update to current displacements
        dw = solution(getDoF(dofm,obj.fldFrac));
        stateCurr.fractureJump(ents) = stateCurr.fractureJump(ents) + dw;

        setState(obj,stateCurr);

        % update traction using penalty approach
        updateTraction(obj);

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
      trac = s.traction(:);

      nCellData = 3;
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      % Displacements
      cellStr(1).name = 'fractureJump';
      cellStr(1).data = [jump(1:3:end), jump(2:3:end), jump(3:3:end)];

      cellStr(2).name = 'traction';
      cellStr(2).data = [trac(1:3:end), trac(2:3:end), trac(3:3:end)];

      cellStr(3).name = 'fractureState';
      as = ContactMode.integer(obj.activeSet.curr);
      cellStr(3).data = as;

      % plot directly into the domain vtm block
      blk = obj.domain.vtmBlock;
      obj.domain.outstate.writeVTKfile(blk,'EmbeddedFractures',obj.fractureMesh,...,
        time,[],[],[],cellStr)
   
    end

   
  end




  methods (Access=private)

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
      initializeGrid(fMesh);
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




    function updateTraction(obj)

      s = getState(obj);
      sOld = getStateOld(obj);
      %sIni = getStateInit(obj);
      jump = s.fractureJump;
      deltaJump = jump - sOld.fractureJump;

      penN = obj.penalty_n;
      penT = obj.penalty_t;
      frac = obj.fractureMesh.surfaces;

      for i = 1:frac.num
        dofW = getLocalDoF(obj.domain.dofm,obj.fldFrac,i);
        t = s.traction(dofW);
        tOld = sOld.traction(dofW);
        j = jump(dofW);
        dj = deltaJump(dofW);
        fId = frac.fracId(i);

        switch obj.activeSet.curr(i)
          case ContactMode.stick
            t(1) = tOld(1) + penN * dj(1);
            t(2:3) = tOld(2:3) + penT * dj(2:3);

          case {ContactMode.slip, ContactMode.newSlip}
            t(1) = tOld(1) + penN * j(1);
            tauLim = obj.cohesion(fId) - t(1)*tan(obj.phi(fId));   % using the updated or not?
            %
            t(2:3) = tauLim * (dj(2:3)/norm(dj(2:3)));
          case ContactMode.open
            t(:) = 0;
        end

        s.traction(dofW) = t;

      end

      setState(obj,s);

      % stick traction might need the initial one...?
      %s.traction = s.traction + sIni.traction;

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



    function dtdg = computeDerTractionGap(obj,i,slip,tN)

      % slip: 2x1 array with tangential gap increment
      state = obj.activeSet.curr(i);
      dtdg = zeros(3,3);
      fId = obj.fractureMesh.surfaces.fracId(i);
      phiF = obj.phi(fId);
      cF = obj.cohesion(fId);
      tauLim = cF - tN*tan(phiF);

      switch state
        case ContactMode.open
          return
        case ContactMode.stick
          dtdg = diag([obj.penalty_n,obj.penalty_t,obj.penalty_t]);
        case {ContactMode.slip,ContactMode.newSlip}
          dtdg(1,1) = obj.penalty_n;
          if norm(slip) > obj.activeSet.tol.sliding
            slipNorm = norm(slip);
            dtdg([2 3],1) = - obj.penalty_n*tan(phiF)*(slip/slipNorm);
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

