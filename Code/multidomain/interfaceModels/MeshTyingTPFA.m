classdef MeshTyingTPFA < InterfaceSolver

  % Finite volume mesh tying between non conforming interfaces

  % One interface pressure per active slave surface. Half transmissibilities
  % are evaluated on intersections; flux continuity is summed per slave face.
  % Now only for SinglePhaseFlow without gravity

  properties
    neigh                     % volume-cell neighbors; columns indexed by MortarSide
    faceArea                  % area of the intersection polygons
    faceCenter                % centroid of the intersection polygons
    halfTrans                 % Nface x 2 half transmissibilities; columns indexed by MortarSide
    flowModels = cell(2,1)    % handle to Flow models in master and slave side
    faceNormals               % normals of the intersection polygons
    intersectionToSlave       % intersection row -> local slave surface / pressure DoF

  end

  methods

    function obj = MeshTyingTPFA(id,domains,inputStruct)
      
      obj@InterfaceSolver(id,domains,inputStruct);

    end

    function registerInterface(obj,input)

      % check validity of the model
      validateInterface(obj);

      % One pressure on each retained slave face, shared by its intersections.
      obj.nMult = obj.grids(MortarSide.slave).surfaces.num;
      obj.multiplierLocation = entityField.surface;
      state = getState(obj);
      state.interfacePressure = zeros(obj.nMult,1);
      setState(obj,state);
      setStateOld(obj,state);
      setStateInit(obj,state);

    end


    function updateState(obj,du)

      p = getState(obj,"interfacePressure");
      setState(obj,p + du(1:obj.nMult),"interfacePressure");

    end


    function assembleConstraint(obj,varargin)

      % set up half face transmissibilities and compute rhs fluxes
      if isempty(obj.halfTrans)
        computeHalfTrans(obj);
      end

      % we are assuming SINGLEPHASEFLOW for now
      % consider replacing with a 
      % flowSolver{side}.computeMobility() method valid for
      % generic flow model
      % mob is an array with size nF x 2
      mob = zeros(2,1);
      mob(1) = 1/obj.domains(1).materials.getFluid().getDynViscosity();
      mob(2) = 1/obj.domains(2).materials.getFluid().getDynViscosity();

      % reset the jacobian blocks
      obj.setJmu(MortarSide.slave, []);
      obj.setJmu(MortarSide.master, []);
      obj.setJum(MortarSide.slave, []);
      obj.setJum(MortarSide.master, []);

      assembleMatrices(obj,mob); 
      [rhsMaster,rhsSlave,rhsInterf] = computeRhs(obj,mob);

      addRhs(obj,MortarSide.master,rhsMaster);
      addRhs(obj,MortarSide.slave,rhsSlave);
      obj.rhsConstraint = rhsInterf;

    end


    function assembleMatrices(obj,mob)

      % assemble half transmissibilities

      Tff = sparse(obj.nMult,obj.nMult);

      for side = [MortarSide.master,MortarSide.slave]

        s = getSide(side);

        hT = mob(s).*obj.halfTrans(:,s);

        dofm = obj.getDoFManager(side);

        idPressure = dofm.getVariableId(obj.coupledVariables);

        nDoF = dofm.getNumbDoF(obj.coupledVariables);

        neighs = getDoFManager(obj,side).getLocalDoF(idPressure,obj.neigh(:,s));

        % Repeated slave indices sum all intersections of the same face.
        Tfu = sparse(obj.intersectionToSlave,neighs,hT,obj.nMult,nDoF);
        Tuu = sparse(neighs,neighs,hT,nDoF,nDoF);

        addJmu(obj,side,-Tfu);
        addJum(obj,side,-Tfu');
        obj.flowModels{s}.domain.J{idPressure,idPressure} = ...
          obj.flowModels{s}.domain.J{idPressure,idPressure} + Tuu;

        Tff = Tff + sparse(obj.intersectionToSlave,...
          obj.intersectionToSlave,hT,obj.nMult,obj.nMult);
        
      end

      obj.Jconstraint = Tff;

    end

    function [rhsM,rhsS,rhsI] = computeRhs(obj,mob)

      pFace = getState(obj,"interfacePressure");
      pI = pFace(obj.intersectionToSlave);
      rhs = cell(2,1);
      rhsI = zeros(obj.nMult,1);

      for side = [MortarSide.slave,MortarSide.master]
        dofm = getDoFManager(obj,side);
        idPressure = dofm.getVariableId(obj.coupledVariables);
        nDoF = dofm.getNumbDoF(obj.coupledVariables);
        cells = obj.neigh(:,side);
        dofs = dofm.getLocalDoF(idPressure,cells);
        p = getState(obj.domains(side),obj.coupledVariables);
        hT = mob(side).*obj.halfTrans(:,side);
        q = hT.*(p(cells)-pI); % outward flux from this volume cell
        rhs{side} = accumarray(dofs,q,[nDoF,1]);
        rhsI = rhsI - accumarray(obj.intersectionToSlave,q,[obj.nMult,1]);
      end

      rhsM = rhs{MortarSide.master};
      rhsS = rhs{MortarSide.slave};

    end


    function writeSolution(obj,fac,tID)

      p = obj.state.interpolate(fac,"interfacePressure");

      obj.outstate.results(tID).interfacePressure = p;

    end


    function [surfaceStr,pointStr] = writeVTK(obj,fac,varargin)

      % The pressure ordering is exactly the retained slave surface ordering.
      surfaceStr = struct('name','interfacePressure',...
        'data',obj.state.interpolate(fac,"interfacePressure"));
      pointStr = [];

    end


    function writeMatFile(obj,fac,tID)

      p = obj.state.interpolate(fac,"interfacePressure");
      obj.outstate.matFile(tID).interfacePressure = p;

    end


  end


  methods (Access = protected)

    function setMortarInterface(obj,params)

      % force segment-based mortar processing
      quad.type = "SegmentBasedQuadrature";
      quad.gaussOrder = 1;
      params.Quadrature = quad;
      params.multiplierType = "P0";

      setMortarInterface@InterfaceSolver(obj,params);

      % Quadrature pairs contain LOCAL surface IDs, one row per subtriangle.
      % Retain intersection geometry, but share pressure on each slave face.
      pairs = obj.quadrature.interfacePairs;
      [surfacePairs,~,polygonId] = unique(pairs,'rows');
      obj.intersectionToSlave = surfacePairs(:,MortarSide.slave);
      obj.multiplierLocation = entityField.surface;
      nConn = size(surfacePairs,1);
      obj.neigh = zeros(nConn,2);
      obj.faceArea = zeros(nConn,1);
      obj.faceCenter = zeros(nConn,3);
      obj.faceNormals = zeros(nConn,3);
      obj.halfTrans = [];

      if nConn == 0
        return
      end

      slave = MortarSide.slave;
      surfSlave = obj.grids(slave).surfaces;
      firstMoment = zeros(size(pairs,1),3);
      triArea = zeros(size(pairs,1),1);

      % Map local Gauss coordinates to physical coordinates BEFORE averaging.
      % The weights already include the subtriangle Jacobian (physical area).
      for vtk = surfSlave.vtkTypes
        elem = getElement(obj,vtk,slave);
        triIds = find(surfSlave.VTKType(pairs(:,slave)) == vtk);
        for k = triIds'
          is = pairs(k,slave);
          xi = obj.quadrature.getSlaveGPCoords(k);
          x = elem.computeBasisF(xi)*elem.getElementCoords(is);
          w = obj.quadrature.getIntegrationWeights(k);
          w = w(:);
          triArea(k) = sum(w);
          firstMoment(k,:) = w'*x;
        end
      end

      obj.faceArea = accumarray(polygonId,triArea,[nConn,1]);
      assert(all(isfinite(obj.faceArea) & obj.faceArea > 0),...
        "MeshTyingTPFA: invalid intersection area");
      for d = 1:3
        obj.faceCenter(:,d) = accumarray(polygonId,firstMoment(:,d),...
          [nConn,1])./obj.faceArea;
      end

      % Surface grids retain faceId through extraction and local renumbering.
      % Convert to parent volume cells, using the enum for column ordering.
      for side = [MortarSide.slave,MortarSide.master]
        faceId = obj.grids(side).surfaces.faceId(surfacePairs(:,side));
        faces = obj.domains(side).grid.faces;
        assert(all(faces.neighbors(faceId,2) == 0),...
          "MeshTyingTPFA requires boundary faces on both sides");
        obj.neigh(:,side) = faces.neighbors(faceId,1);
      end

      % The shared normal points outward from the slave volume.
      obj.faceNormals = surfSlave.normal(surfacePairs(:,slave),:);

    end


    function computeHalfTrans(obj)

      % half face transmissibility 
      % the property is a Nfaces x 2 matrix with the half face
      % transmissibility of each neigboring cell
      r = [1, 1, 1, 2, 2, 2, 3, 3, 3];
      c = [1, 2, 3, 1, 2, 3, 1, 2, 3];

      obj.halfTrans = zeros(length(obj.faceArea),2);

      
      % sign to ensure that the scalar product between c and normal is positive
      % remember: normal points outward the slave domain
      sgn = zeros(1,2);
      sgn(MortarSide.master) = -1;
      sgn(MortarSide.slave) = +1;

      for side = [MortarSide.master,MortarSide.slave]

        hT = zeros(size(obj.halfTrans,1),1);

        s = getSide(side);
        mesh = obj.domains(s).grid;
        L = obj.faceCenter - mesh.cells.center(obj.neigh(:,s),:);
        L = sgn(s)*L;
        KMat = zeros(mesh.cells.nTag,9);
        for i=1:mesh.cells.nTag
          KMat(i,:) = obj.domains(s).materials.getMaterial(i).PorousRock.getPermVector();
        end

        for k=1:length(r)
          hT = hT + L(:,r(k)) .* KMat(mesh.cells.tag(obj.neigh(:,s)),k) .* obj.faceNormals(:,c(k));
        end

        hT = obj.faceArea.*hT;
        obj.halfTrans(:,s) = hT./sum(L.*L,2);

      end

    end 




  end

  methods (Access=private)

    function  validateInterface(obj)

      % validate the model and populate the obj.flowModels pointers

      var = obj.coupledVariables;
      assert(strcmp(obj.coupledVariables,"pressure"),...
        "MeshTyingTPFA can only couple pressure field with TPFA");

      for side = [MortarSide.master,MortarSide.slave]

        s = getSide(side);
        fieldLocation = getDoFManager(obj,side).getFieldLocation(var);
        assert(fieldLocation == entityField.cell,...
          "MeshTyingTPFA works only with cell pressure variables in master and slave domains")

        solvNames = obj.domains(s).solverNames;

        gamma = obj.domains(s).materials.getFluid().getSpecificWeight();
        assert(gamma==0.0,"Gravity contribution not yet supported in MeshTyingTPFA")

        obj.flowModels{s} = getPhysicsSolver(obj.domains(s),solvNames); 

      end

    end
  end

  methods (Static)

    function var = getCoupledVariables()
      var = [];
    end

  end
end



