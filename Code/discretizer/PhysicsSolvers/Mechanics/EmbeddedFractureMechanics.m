classdef EmbeddedFractureMechanics < PhysicsSolver

  % solver for embedded tractions implementing the EFEM(0) formulation
  % Cusini et al (2021).

  properties

    cutCells                % list of global index of cells intercepted by fracture
    cutAreas    (:,1)       % the area of the cut cells
    cutCenters  (:,3)  
    cutNormals  (:,3)       % the normal of the cutting plane for the fracture
    cutTang1    (:,3)
    cutTang2    (:,3)
    nCutCells
    contactState = struct("curr",[],"old",[])          
    penalty_n               % penalty parameter for normal direction
    penalty_t               % penalty parameter for tangential direction
    phi                     % friction angle in radians for each fracture
    cohesion                % the cohesion of each fracture
    cutCellToFracture       % map each cut cell to its fracture id
    fractureMesh            % a 2D mesh object with cut cell topology
    outFracture
    areaTol = 1e-6;         % minimum area of a fracture element

  end

  properties (Access = private)
    fieldId
  end

  methods (Access = public)

    function obj = EmbeddedFractureMechanics(domain)

      % call physicsSolver constructor
      obj@PhysicsSolver(domain);

    end

    function registerSolver(obj,solverInput)

      obj.penalty_n = getXMLData(solverInput,[],"penaltyNormal");
      obj.penalty_t = getXMLData(solverInput,[],"penaltyTangential");

      defineFractures(obj,solverInput);

      % register nodal displacements on target regions
      obj.dofm.registerVariable(obj.getField(),entityField.cell,3,"nEntities",obj.nCutCells);

      % store the id of the field in the degree of freedom manager
      obj.fieldId = obj.dofm.getVariableId(obj.getField());

      % initialize the state object
      initState(obj);

      obj.outFracture = copy(obj.domain.outstate);

      outFileName = getXMLData(solverInput,"fractureOutput","outputFile");
      obj.outFracture.VTK = VTKOutput(obj.fractureMesh,outFileName);

    end

    function assembleSystem(obj,dt)
      % compute the displacements matrices and rhs in the domain
      fldMech = obj.dofm.getVariableId(Poromechanics.getField());
      [Juw,Jwu,Jww,rhsW] = computeJacobianAndRhs(obj,dt);

      obj.domain.J{obj.fieldId,fldMech} = Jwu;
      obj.domain.J{fldMech,obj.fieldId} = Juw;
      obj.domain.J{obj.fieldId,obj.fieldId} = Jww;
      obj.domain.rhs{obj.fieldId} = rhsW;

    end


    function [Kuw,Kwu,Kww,rhsW] = computeJacobianAndRhs(obj,dt)

      n1 = sum((obj.mesh.nDim^2)*(obj.mesh.cellNumVerts(obj.cutCells)*obj.nCutCells));
      n2 = sum((obj.mesh.nDim^2)*obj.nCutCells^2);
      nDofU = obj.dofm.getNumbDoF(Poromechanics.getField()); 
      nDofW = obj.dofm.getNumbDoF(obj.fieldId); 

      asbKuw = assembler(n1,nDofU,nDofW);
      asbKwu = assembler(n1,nDofW,nDofU);
      asbKww = assembler(n2,nDofW,nDofW);

      rhsW = zeros(nDofW,1);

      % get the state object
      s = getState(obj);
      sOld = getStateOld(obj);

      cell2stress = getPhysicsSolver(obj.domain,"Poromechanics").cell2stress;

      fldMech = obj.dofm.getVariableId(Poromechanics.getField());

      for i = 1:obj.nCutCells

        % compute local terms and assemble to requested matrices
        cellId = obj.cutCells(i);

        l = cell2stress(cellId);

        % compute B matrix
        vtkId = obj.mesh.cellVTKType(cellId);
        elem = getElement(obj.elements,vtkId);
        nG = elem.GaussPts.nNode;
        [N,dJWeighed] = getDerBasisFAndDet(elem,cellId,1);
        B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
        B(elem.indB(:,2)) = N(elem.indB(:,1));

        % compute Bw matrix (compatibility operator, 6x3)
        Bw = computeCompatibilityMatrix(obj,i,N);

        % compute E matrix (equilibrium operator, 6x3)
        E = computeEquilibriumOperator(obj,i);

        % compute constituvie tensor
        [D, sigma, ~] = obj.materials.updateMaterial( ...
          obj.mesh.cellTag(cellId), ...
          sOld.data.stress(l+1:l+nG,:), ...
          s.data.strain(l+1:l+nG,:), ...
          dt, sOld.data.status(l+1:l+nG,:), cellId, s.t);

        KuwLoc = Poromechanics.computeKloc(B,D,Bw,dJWeighed);
        KwuLoc = Poromechanics.computeKloc(E,D,B,dJWeighed);
        KwwLoc = Poromechanics.computeKloc(E,D,Bw,dJWeighed);

        dtdg = computeDerTractionGap(obj,i);

        KwwLoc = KwwLoc - dtdg*obj.cutAreas(i);

        % grab degree of freedom
        nodes = obj.mesh.cells(cellId,:);
        uDof = obj.dofm.getLocalDoF(fldMech,nodes);
        wDof = obj.dofm.getLocalDoF(obj.fieldId,i);

        % assemble local contributions
        asbKuw.localAssembly(uDof,wDof,KuwLoc);
        asbKwu.localAssembly(wDof,uDof,-KwuLoc);
        asbKww.localAssembly(wDof,wDof,-KwwLoc);

        % assemble rhsW (use computed stress tensor)
        sigma = reshape(sigma',6,1,nG);
        trac = s.data.traction(wDof);
        r1 = trac*obj.cutAreas(i);
        fTmp = pagemtimes(E,'ctranspose',sigma,'none');
        fTmp = fTmp.*reshape(dJWeighed,1,1,[]);
        r2 = sum(fTmp,3);
        rhsLoc = - r1 - r2;
        rhsW(wDof) = rhsW(wDof) + rhsLoc; 

      end

      Kuw = asbKuw.sparseAssembly();
      Kwu = asbKwu.sparseAssembly();
      Kww = asbKww.sparseAssembly();

    end


    function initState(obj)
      % add poromechanics fields to state structure
      state = getState(obj);
      state.data.(obj.getField()) = zeros(3*obj.nCutCells,1);
      state.data.traction = zeros(3*obj.nCutCells,1);
      obj.contactState.curr = repmat(ContactMode.stick,obj.nCutCells,1);
      obj.contactState.old = obj.contactState.curr;

    end

    function advanceState(obj)
      % Set converged state to current state after newton convergence
      obj.contactState.old = obj.contactState.curr;
    end

    function updateState(obj,solution)


      % Update state structure with last solution increment
      ents = obj.dofm.getActiveEntities(obj.fieldId,1);

      stateCurr = obj.getState();
      stateOld = obj.getStateOld();



      if nargin > 1
        % apply newton update to current displacements
        dw = solution(getDoF(obj.dofm,obj.fieldId));
        stateCurr.data.(obj.getField())(ents) = stateCurr.data.(obj.getField())(ents) + ...
          dw;

        % update traction using penalty approach
        updateTraction(obj);


      end

      cell2stress = getPhysicsSolver(obj.domain,"Poromechanics").cell2stress;
       
      % jump increment at current iteration
      w = stateCurr.data.(obj.getField()) - stateOld.data.(obj.getField()); 

      % Enhance straint with fracture contribution
      for i = 1:obj.nCutCells
        el = obj.cutCells(i);
        l = cell2stress(el);
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        nG = elem.GaussPts.nNode;
        N = getDerBasisFAndDet(elem,el,2);
        Bw = computeCompatibilityMatrix(obj,i,N);
        jump = w(getLocalDoF(obj.dofm,obj.fieldId,i));
        stateCurr.data.strain(l+1:l+nG,:) = stateCurr.data.strain(l+1:l+nG,:) + ...
          reshape(pagemtimes(Bw,jump),6,nG)';
      end
    end

    function [avStress,avStrain] = finalizeState(obj,stateIn)
      % compute cell average values of stress and strain
      avStress = zeros(obj.mesh.nCells,6);
      avStrain = zeros(obj.mesh.nCells,6);

      l = 0;
      for el = 1:obj.mesh.nCells
        dof = getDoFID(obj.mesh,el);
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        nG = elem.GaussPts.nNode;
        vol = obj.mesh.cellVolume(el);
        [N,dJWeighed] = getDerBasisFAndDet(elem,el,1);
        B = zeros(6,elem.nNode*obj.mesh.nDim,nG);
        B(elem.indB(:,2)) = N(elem.indB(:,1));
        avStress(el,:) = sum(diag(dJWeighed)* ...
          stateIn.data.stress((l+1):(l+nG),:))/vol;
        dStrain = pagemtimes(B,stateIn.data.displacements(dof));
        dStrain = dStrain.*reshape(dJWeighed,1,1,[]);
        avStrain(el,:) = sum(dStrain,3)/vol;
        l = l + nG;
      end
    end

    function applyBC(obj,bcId,t)

      if ~BCapplies(obj,bcId)
        return
      end

      % get bcDofs and bcVals
      [bcDofs,bcVals] = getBC(obj,bcId,t);

      bcType = obj.bcs.getType(bcId);

      switch bcType
        case 'Dirichlet'
          applyDirBC(obj,bcId,bcDofs);
        case {'Neumann','VolumeForce'}
          applyNeuBC(obj,bcId,bcDofs,bcVals);
        otherwise
          error("Error in %s: Boundary condition type '%s' is not " + ...
            "available in %s",class(obj),bcType);
      end
      
    end


    function [dof,vals] = getBC(obj,bcId,t)
      %
      dof = obj.getBCdofs(bcId);
      vals = obj.getBCVals(bcId,t);
    end


    function applyDirVal(obj,bcId,t)

      bcVar = obj.bcs.getVariable(bcId);

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

    function writeMatFile(obj,fac,tID)

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
      cellData(1).data = double(reshape(ismember(1:obj.mesh.nCells,obj.cutCells),[],1));
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

      nCellData = 2;
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      % Displacements
      cellStr(1).name = 'fractureJump';
      cellStr(1).data = [jump(1:3:end), jump(2:3:end), jump(3:3:end)];

      cellStr(2).name = 'traction';
      cellStr(2).data = [trac(1:3:end), trac(2:3:end), trac(3:3:end)];

      % plot with the fracture output object
      obj.outFracture.VTK.writeVTKFile(time, [], [], [], cellStr);

    end

    function finalizeOutput(obj)
      obj.outFracture.finalize();
    end

  end

  methods (Access=private)

    function defineFractures(obj,input)

      fractureStruct = input.Fracture;

      nFractures = numel(fractureStruct);

      [edgesTopol,cellToEdges] = getEdgeTopology(obj);

      nE = size(edgesTopol,1);

      tol = 1e-8;

      countCell = 0;
      obj.cutCells = zeros(obj.mesh.nCells*nFractures,1);

      obj.cohesion = zeros(nFractures,1);
      obj.phi = zeros(nFractures,1);


      obj.fractureMesh = Mesh();
      msh = obj.fractureMesh;

      nV = 0;


      for f = 1:nFractures

        % read the fracture geometry
        normal = getXMLData(fractureStruct(f),[],'normal');
        origin = getXMLData(fractureStruct(f),[],'origin');
        dims = getXMLData(fractureStruct(f),[],'dimensions');
        lVec = getXMLData(fractureStruct(f),[],'lengthVec');
        wVec = getXMLData(fractureStruct(f),[],'widthVec');

        obj.cohesion(f) = getXMLData(fractureStruct(f),[],'cohesion');
        obj.phi(f) = getXMLData(fractureStruct(f),[],'frictionAngle');


        assert(all([abs(lVec*normal') < 1e-8,...
              abs(wVec*normal')<1e-8,...
              abs(lVec*wVec')<1e-8]),...
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
        distVec = obj.mesh.coordinates - origin;
        nVec = reshape(normal/norm(normal),[],1);
        R = InterfaceMesh.computeRot(nVec);
        tVec1 = R(:,2);
        tVec2 = R(:,3);

        sign = distVec * nVec > 0;

        % check edges with nodes in different locations
        isEdgeCrossed = sum(sign(edgesTopol),2) == 1;
       
        % compute intersection points of crossed edges
        x0 = obj.mesh.coordinates(edgesTopol(isEdgeCrossed,1),:);
        x1 =  obj.mesh.coordinates(edgesTopol(isEdgeCrossed,2),:);
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

        obj.cutCells(countCell+1:countCell+nC) = newCutCells;

        obj.cutNormals(countCell+1:countCell+nC,:) = repmat(nVec',nC,1);
        obj.cutTang1(countCell+1:countCell+nC,:) = repmat(tVec1',nC,1);
        obj.cutTang2(countCell+1:countCell+nC,:) = repmat(tVec2',nC,1);
        obj.cutCellToFracture(countCell+1:countCell+nC,:) = f;

        % preallocate number of nodes
        surfs = zeros(nC,6);
        surfNumVerts = zeros(nC,1);

        % loop over cut cell and compute geometry
        for ic = 1:nC
          cellEdges = cellToEdges(newCutCells(ic),:);
          isEdgeCut = logical(m(newCutCells(ic),:));
          cutEdges =  cellEdges(isEdgeCut);
          cutCellVertices = intersections(cutEdges,:);

          idx = orderPointsCCW(cutCellVertices,normal);

          surfs(ic,1:numel(cutEdges)) = cutEdges(idx);
          surfNumVerts(ic) = numel(cutEdges);
              
          [obj.cutCenters(countCell+ic,:),obj.cutAreas(countCell+ic)] = ...
            computePolygonGeometry(cutCellVertices,nVec');

        end

        countCell = countCell + nC;

        % finalize the mesh for the current fracture
        [~,~,id] = unique(surfs(:));

        surfs = nV + reshape(id,[],6);
        surfs = surfs - 1;

        nV = sum(id > 1);

        msh.surfaces = [msh.surfaces; surfs];
        msh.surfaceNumVerts = [msh.surfaceNumVerts; surfNumVerts];
        msh.coordinates = [msh.coordinates; xInt(isInPlane,:)];



      end

        obj.cutCells = obj.cutCells(1:countCell);
        obj.cutNormals = obj.cutNormals(1:countCell,:);
        obj.cutCenters = obj.cutCenters(1:countCell,:);
        obj.cutAreas = obj.cutAreas(1:countCell);

        % finalize fracture mesh

        % discard too small fractures
        id = obj.cutAreas > obj.areaTol;
        obj.nCutCells = sum(id);
        obj.cutCells = obj.cutCells(id);
        obj.cutCenters = obj.cutCenters(id,:);
        obj.cutAreas = obj.cutAreas(id);
        obj.cutNormals = obj.cutNormals(id,:);
        obj.cutTang1 = obj.cutTang1(id,:);
        obj.cutTang2 = obj.cutTang2(id,:);

        surfs = msh.surfaces(id,:);
      
        [u,~,id2] = unique(surfs(:));

        surfs = reshape(id2,[],6);
        msh.surfaces = surfs - 1;

        msh.surfaceNumVerts = msh.surfaceNumVerts(id);
        
        msh.coordinates = msh.coordinates(u(2:end),:);

        % coordRound = round(msh.coordinates/1e-7);
        % [nodesUnique, ia, ic] = unique(coordRound, 'rows', 'stable')

        msh.nNodes = size(msh.coordinates,1);
        msh.nSurfaces = size(msh.surfaces,1);
        msh.surfaceVTKType = 7*ones(msh.nSurfaces,1);

    end


    function [edges,c2e] = getEdgeTopology(obj)

      assert(isscalar(unique(obj.mesh.cellVTKType)),...
        "EFEM implemented only for mesh with uniform element shapes");


      if obj.mesh.cellVTKType(1) == 12

        eLoc = [ ...
          1 2; 2 3; 3 4; 4 1; ...
          5 6; 6 7; 7 8; 8 5; ...
          1 5; 2 6; 3 7; 4 8 ];

        nCells = obj.mesh.nCells;
        nEdgesLoc = size(eLoc,1);

        cells = obj.mesh.cells;

        % Node indices of all edges (stacked)
        e1 = cells(:, eLoc(:,1));   % (nCells x 12)
        e2 = cells(:, eLoc(:,2));   % (nCells x 12)

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

      % check

      s = getState(obj);
      sOld = getStateOld(obj);
      jump = s.data.(obj.getField());
      deltaJump = jump - sOld.data.(obj.getField()); 

      penN = obj.penalty_n;
      penT = obj.penalty_t;

      for i = 1:obj.nCutCells
        dofW = getLocalDoF(obj.dofm,obj.fieldId,i);
        j = jump(dofW);
        dj = deltaJump(dofW);

        switch obj.contactState.curr(i)
          
          case ContactMode.stick
            s.data.traction(dofW(1)) = ...
              s.data.traction(dofW(1)) + penN * j(1);
            s.data.traction(dofW(2:3)) = ...
              s.data.traction(dofW(2:3)) + penT * j(2:3);

          case ContactMode.slip
             s.data.traction(dofW(1)) = ...
               s.data.traction(dofW(1)) + penN * j(1);
             tN = s.data.traction(dofW(1));
            tauLim = obj.cohesion - tN*tan(obj.phi);
            %
            s.data.traction(dofW(2:3)) = tauLim * (dj(2:3)/norm(dj(2:3)));
            
          case ContactMode.open
            s.data.traction(dofW(:)) = 0;

        end


      end
    end


    function H = computeHeaviside(obj,i)

      cellId = obj.cutCells(i);
      coords = obj.mesh.coordinates(obj.mesh.cells(cellId,:),:);
      dist = obj.cutCenters(i) - coords;
      n = obj.cutNormals(i,:)';
      H = double(dist*n > 0);

    end

    function Bw = computeCompatibilityMatrix(obj,i,N)

      H = computeHeaviside(obj,i);
      v = - sum(N.*H',2);
      n = obj.cutNormals(i,:);
      m1 = obj.cutTang1(i,:);
      m2 = obj.cutTang2(i,:);

      v = permute(v,[2 1 3]);


      sym_n_dyad_v = obj.sym_AiBj_plus_AjBi(n,v);
      sym_m1_dyad_v = obj.sym_AiBj_plus_AjBi(m1,v);
      sym_m2_dyad_v = obj.sym_AiBj_plus_AjBi(m2,v);

      Bw = [sym_n_dyad_v, sym_m1_dyad_v, sym_m2_dyad_v];
    end

    function E = computeEquilibriumOperator(obj,i)

      n = obj.cutNormals(i,:);
      m1 = obj.cutTang1(i,:)';
      m2 = obj.cutTang2(i,:)';

      sym_n_dyad_n = obj.sym_AiBj_plus_AjBi(n,n);
      sym_m1_dyad_n = obj.sym_AiBj_plus_AjBi(m1,n);
      sym_m2_dyad_n = obj.sym_AiBj_plus_AjBi(m2,n);

      A = obj.cutAreas(i);
      V = obj.mesh.cellVolume(obj.cutCells(i));

      E = - (A/V) * [sym_n_dyad_n, sym_m1_dyad_n, sym_m2_dyad_n];

    end

    function dtdg = computeDerTractionGap(obj,i,slip)

      % slip: 2x1 array with tangential gap increment
      state = obj.contactState.curr(i);
      dtdg = zeros(3,3);
      tN = obj.domain.state.data.traction(3*i-2);
      tauLim = obj.cohesion - tN*tan(obj.phi);

      switch state
        case ContactMode.open
          return
        case ContactMode.stick
          dtdg = -diag([obj.penalty_n,obj.penalty_t,obj.penalty_t]);
        case ContactMode.slip
          dtdg(1,1) = obj.penalty_n;
          slipNorm = norm(slip);
          dtdg([2 3],1) = - obj.penalty_n*tan(obj.phi)*(slip/slipNorm);
          dtdg([2 3],[2 3]) = tauLim * ...
            (slipNorm^2*eye(3) - slip * slip')/slipNorm^3; 

      end
     
    end

  end

  methods (Static)
    
    function [cellStr,pointStr] = buildPrintStruct(disp,stress,strain)

      nCellData = 2;
      nPointData = 1;
      pointStr = repmat(struct('name', 1, 'data', 1), nPointData, 1);
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      % Displacements
      pointStr(1).name = 'displacements';
      pointStr(1).data = [disp(1:3:end) disp(2:3:end) disp(3:3:end)];
      % pointStr(1).name = 'ux';
      % pointStr(1).data = disp(1:3:end);
      % pointStr(2).name = 'uy';
      % pointStr(2).data = disp(2:3:end);
      % pointStr(3).name = 'uz';
      % pointStr(3).data = disp(3:3:end);
      %

      % Stress
      cellStr(1).name = 'stress';
      cellStr(1).data = stress;
      cellStr(2).name = 'strain';
      cellStr(2).data = strain;
      % cellStr(3).name = 'sz';
      % cellStr(3).data = stress(:,3);
      % cellStr(4).name = 'txy';
      % cellStr(4).data = stress(:,4);
      % cellStr(5).name = 'tyz';
      % cellStr(5).data = stress(:,5);
      % cellStr(6).name = 'txz';
      % cellStr(6).data = stress(:,6);
      % %
      % % Strain
      % cellStr(7).name = 'ex';
      % cellStr(7).data = strain(:,1);
      % cellStr(8).name = 'ey';
      % cellStr(8).data = strain(:,2);
      % cellStr(9).name = 'ez';
      % cellStr(9).data = strain(:,3);
      % cellStr(10).name = 'gxy';
      % cellStr(10).data = strain(:,4);
      % cellStr(11).name = 'gyz';
      % cellStr(11).data = strain(:,5);
      % cellStr(12).name = 'gxz';
      % cellStr(12).data = strain(:,6);
    end


    function vSym = sym_AiBj_plus_AjBi(a,b)
      % compute symmetric dyadic product of two 3x3 tensor 
      % return result into a 6x1 voigt array

      a = reshape(a,3,1,[]);
      b = reshape(b,3,1,[]);

      v = pagemtimes(a,'none',b,'ctranspose') + pagemtimes(b,'none',a,'ctranspose');
      vSym = reshape(v,9,1,[]);
      vSym = vSym([1 5 9 4 8 7],:,:);
      vSym(1:3,:,:) = 0.5*vSym(1:3,:,:);

    end

    function out = getField()
      out = "fractureJump";
    end

end

end

