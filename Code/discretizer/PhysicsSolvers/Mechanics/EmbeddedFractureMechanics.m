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

  end

  properties (Access = private)
    fieldId
  end

  methods (Access = public)

    function obj = EmbeddedFractureMechanics(domain)

      % call physicsSolver constructor
      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,varargin)

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
      dofm.registerVariable(obj.getField(),entityField.cell,3,"nEntities",nCutCells);

      % store the id of the field in the degree of freedom manager
      obj.fieldId = dofm.getVariableId(obj.getField());

      % initialize the state object
      initState(obj);

      initializeActiveSet(obj,nCutCells,params.ActiveSet);

    end

    function assembleSystem(obj,dt)
      % compute the displacements matrices and rhs in the domain

      fldMech = obj.domain.dofm.getVariableId(Poromechanics.getField());
      [Juw,Jwu,Jww,rhsW] = computeJacobianAndRhs(obj,dt);

      obj.domain.J{obj.fieldId,fldMech} = Jwu;
      obj.domain.J{fldMech,obj.fieldId} = Juw;
      obj.domain.J{obj.fieldId,obj.fieldId} = Jww;
      obj.domain.rhs{obj.fieldId} = rhsW;

    end


    function [Kuw,Kwu,Kww,rhsW] = computeJacobianAndRhs(obj,dt)

      dofm = obj.domain.dofm;
      f = obj.fractureMesh.surfaces;
      cells = obj.grid.cells;
      n1 = sum((obj.grid.nDim^2)*(cells.numVerts(f.cutCells)*f.num));
      n2 = sum((obj.grid.nDim^2)*f.num^2);
      nDofU = dofm.getNumbDoF(Poromechanics.getField());
      nDofW = dofm.getNumbDoF(obj.fieldId);

      asbKuw = assembler(n1,nDofU,nDofW);
      asbKwu = assembler(n1,nDofW,nDofU);
      asbKww = assembler(n2,nDofW,nDofW);

      rhsW = zeros(nDofW,1);

      % get the state object
      s = getState(obj);
      sOld = getStateOld(obj);

      jump = s.data.(obj.getField());

      coordinates = obj.grid.coordinates;
      mech = getPhysicsSolver(obj.domain,"Poromechanics");
      cell2stress = mech.cell2stress;

      fldMech = dofm.getVariableId(Poromechanics.getField());

      % only hexa for now
      elem = Hexahedron(obj.grid,'gaussOrder',mech.getGaussOrder);
      nG = getNumbGaussPts(elem);

      topol = obj.grid.getCellNodes(f.cutCells);

      for i = 1:f.num

        % compute local terms and assemble to requested matrices
        cellId = f.cutCells(i);

        l = cell2stress(cellId);

        nodes = topol(i,:);
        coord = coordinates(nodes,:);

        % compute B matrix
        %dof = dofId(nodes,3);

        [gradN,dJw] = getDerBasisFAndDet(elem,coord);
        B = elem.getStrainMatrix(gradN);

        % compute Bw matrix (compatibility operator, 6x3)
        Bw = computeCompatibilityMatrix(obj,f,i,coord,gradN);

        % compute E matrix (equilibrium operator, 6x3)
        E = computeEquilibriumOperator(obj,f,i);

        % compute constituvie tensor
        [D, sigma, ~] = obj.domain.materials.updateMaterial( ...
          cells.tag(cellId), ...
          sOld.data.stress(l:(l+nG-1),:), ...
          s.data.strain(l:(l+nG-1),:), ...
          dt, sOld.data.status(l:(l+nG-1),:), cellId, s.t);

        KuwLoc = Poromechanics.computeKloc(B,D,Bw,dJw);
        KwuLoc = Poromechanics.computeKloc(E,D,B,dJw);
        KwwLoc = Poromechanics.computeKloc(E,D,Bw,dJw);

        % grab degree of freedom
        uDof = dofm.getLocalDoF(fldMech,nodes);
        wDof = dofm.getLocalDoF(obj.fieldId,i);

        dtdg = computeDerTractionGap(obj,i,jump(wDof([2;3])));

        KwwLoc = KwwLoc - dtdg*f.area(i);

        % assemble local contributions
        asbKuw.localAssembly(uDof,wDof,KuwLoc);
        asbKwu.localAssembly(wDof,uDof,KwuLoc);
        asbKww.localAssembly(wDof,wDof,KwwLoc);

        % assemble rhsW (use computed stress tensor)
        sigma = reshape(sigma',6,1,nG);
        trac = s.data.traction(wDof);
        rT = trac*f.area(i);
        fTmp = pagemtimes(E,'ctranspose',sigma,'none');
        fTmp = fTmp.*reshape(dJw,1,1,[]);
        rSigma = sum(fTmp,3);
        rBC = obj.bcTraction(wDof)*f.area(i);
        rhsLoc = rSigma - rT - rBC;
        rhsW(wDof) = rhsW(wDof) + rhsLoc;


      end

      Kuw = asbKuw.sparseAssembly();
      Kwu = asbKwu.sparseAssembly();
      Kww = asbKww.sparseAssembly();

    end


    function initState(obj)
      % add poromechanics fields to state structure
      nCutCells = obj.fractureMesh.surfaces.num;
      state = getState(obj);
      state.data.(obj.getField()) = zeros(3*nCutCells,1);
      state.data.traction = zeros(3*nCutCells,1);
      obj.bcTraction = zeros(3*nCutCells,1);
      obj.activeSet.curr = repmat(ContactMode.open,nCutCells,1);
      obj.activeSet.prev = obj.activeSet.curr;

    end


    function hasConfigurationChanged = updateConfiguration(obj)

      oldActiveSet = obj.activeSet.curr;

      traction = getState(obj,"traction");
      displacementJump = getState(obj,obj.getField());

      f = obj.fractureMesh.surfaces;


      for i = 1:numel(obj.activeSet.curr)

        state = obj.activeSet.curr(i);

        id = DoFManager.dofExpand(i,3);

        t = traction(id);
        g_n = displacementJump(id(1));

        limitTraction = abs(obj.cohesion - tan(obj.phi)*t(1));

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

    function advanceState(obj)
      % Set converged state to current state after newton convergence
      obj.activeSet.old = obj.activeSet.curr;
    end

    function updateState(obj,solution)

      % Update state structure with last solution increment
      dofm = obj.domain.dofm;
      ents = dofm.getActiveEntities(obj.fieldId,1);

      stateCurr = obj.getState();
      stateOld = obj.getStateOld();

      if nargin > 1
        % apply newton update to current displacements
        dw = solution(getDoF(dofm,obj.fieldId));
        stateCurr.data.(obj.getField())(ents) = stateCurr.data.(obj.getField())(ents) + ...
          dw;

        % update traction using penalty approach
        updateTraction(obj);


      end

      mech = getPhysicsSolver(obj.domain,"Poromechanics");
      cell2stress = mech.cell2stress;
      gOrd = mech.getGaussOrder;

      % jump increment at current iteration
      w = stateCurr.data.(obj.getField()) - stateOld.data.(obj.getField());

      f = obj.fractureMesh.surfaces;
      elem = Hexahedron(obj.grid,'gaussOrder',gOrd);
      nG = elem.getNumbGaussPts;

      topol = obj.grid.getCellNodes(f.cutCells);
      coords = obj.grid.coordinates;

      % Enhance strain with fracture contribution
      for i = 1:f.num
        el = f.cutCells(i);
        l = cell2stress(el);
        nodes = topol(i,:);
        coord = coords(nodes,:);
        gradN = getDerBasisFAndDet(elem,coord);
        Bw = computeCompatibilityMatrix(obj,f,i,coord,gradN);
        jump = w(getLocalDoF(dofm,obj.fieldId,i));
        stateCurr.data.strain(l:(l+nG-1),:) = stateCurr.data.strain(l:(l+nG-1),:) + ...
          reshape(pagemtimes(Bw,jump),6,nG)';
      end
    end


    function applyBC(obj,bcId,t)

      if ~BCapplies(obj,bcId)
        return
      end

    end

    function applyDirVal(obj,bcId,t)

      bcVar = obj.domain.bcs.getVariable(bcId);

      if ~strcmp(bcVar,obj.getField())
        return
      end

      [bcDofs,bcVals] = getBC(obj,bcId,t);

      obj.getState().data.displacements(bcDofs) = bcVals;
    end


    function out = isLinear(obj)
      out = false;
    end

    function out = isSymmetric(obj)

      % if the problem is linear, then Poromechanics is symmetric
      out = false;

    end

    function writeSolution(obj,fac,tID)

      jumpOld = getStateOld(obj,obj.getField());
      jumpCurr = getState(obj,obj.getField());

      tractionOld = getStateOld(obj,"traction");
      tractionCurr = getState(obj,"traction");

      obj.domain.outstate.results(tID).(obj.getField()) = jumpCurr*fac+jumpOld*(1-fac);
      obj.domain.outstate.results(tID).traction = tractionCurr*fac+tractionOld*(1-fac);

    end

    function [cellData,pointData] = writeVTK(obj,fac,time)

      cellData = repmat(struct('name', 1, 'data', 1), 1, 1);
      cellData(1).name = 'isCellFractured';
      isCutCell = false(obj.grid.cells.num,1);
      isCutCell(obj.fractureMesh.surfaces.cutCells) = true;
      cellData(1).data = double(isCutCell);
      pointData = [];

      % this method do not return outputs for the 3D mesh grid. instead, it
      % works on a separate 2D polygonal mesh

      s = getState(obj);
      sOld = getStateOld(obj);

      tracCurr = s.data.traction;
      tracOld = sOld.data.traction;

      jumpCurr = s.data.(obj.getField());
      jumpOld = sOld.data.(obj.getField());

      trac = tracCurr*fac + tracOld*(1-fac);
      jump = jumpCurr*fac + jumpOld*(1-fac);

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

    function finalizeOutput(obj)
      obj.outFracture.finalize();
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
      [f.tang1,f.tang2,f.cutCells,f.nFracCells] = deal([]);
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
        x1 =  coords(edgesTopol(isEdgeCrossed,2),:);
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
          cutEdges =  cellEdges(isEdgeCut);
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
        f.nFracCells  = [f.nFracCells; nC];

        % finalize the mesh for the current fracture
        surfs = surfs(1:c);
        [~,~,id] = unique(surfs);
        surfs = nV + id;
        nV = sum(id > 1);

        f.connectivity = [f.connectivity; ArrayOfArrays(surfs,cutNumVerts)];
        f.numVerts = [f.numVerts; cutNumVerts];
        fMesh.coordinates = [fMesh.coordinates; xInt(isInPlane,:)];

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
      jump = s.data.(obj.getField());
      deltaJump = jump - sOld.data.(obj.getField());

      penN = obj.penalty_n;
      penT = obj.penalty_t;

      for i = 1:obj.fractureMesh.surfaces.num
        dofW = getLocalDoF(obj.domain.dofm,obj.fieldId,i);
        j = jump(dofW);
        dj = deltaJump(dofW);

        switch obj.activeSet.curr(i)

          case ContactMode.stick
            s.data.traction(dofW(1)) = penN * j(1);
            s.data.traction(dofW(2:3)) = penT * j(2:3);

          case {ContactMode.slip, ContactMode.newSlip}
            s.data.traction(dofW(1)) = penN * j(1);
            tN = s.data.traction(dofW(1));
            tauLim = obj.cohesion - tN*tan(obj.phi);
            %
            s.data.traction(dofW(2:3)) = tauLim * (dj(2:3)/norm(dj(2:3)));

          case ContactMode.open
            s.data.traction(dofW(:)) = 0;

        end

        %s.data.traction = s.data.traction + s.data.iniTraction;


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

    function dtdg = computeDerTractionGap(obj,i,slip)

      % slip: 2x1 array with tangential gap increment
      state = obj.activeSet.curr(i);
      dtdg = zeros(3,3);
      tN = obj.domain.state.data.traction(3*i-2);
      tauLim = obj.cohesion - tN*tan(obj.phi);

      switch state
        case ContactMode.open
          return
        case ContactMode.stick
          dtdg = diag([obj.penalty_n,obj.penalty_t,obj.penalty_t]);
        case {ContactMode.slip,ContactMode.newSlip}
          dtdg(1,1) = obj.penalty_n;
          if norm(slip) > obj.activeSet.tol.sliding
            slipNorm = norm(slip);
            dtdg([2 3],1) = - obj.penalty_n*tan(obj.phi)*(slip/slipNorm);
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
      out = "fractureJump";
    end

  end

end

