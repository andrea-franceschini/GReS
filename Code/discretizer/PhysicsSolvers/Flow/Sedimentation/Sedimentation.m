classdef Sedimentation < PhysicsSolver
  % Sedimentation model subclass

  properties
    grid gridForSedimentation

    sedimentHistory = SedimentsMap()
    sedimentAcc (:,:)
    heighControl double = 0.1

    facesNeigh (:,2) uint64
    trans (:,1)

    H (:,:)
    P (:,:)
    rhsGrav (:,1)         % gravity contribution to rhs
  end

  properties (Access = protected)
    fieldId
  end

  properties (Access =private)
    grow logical = false
    nmat uint16
  end

  methods (Access = public)
    function obj = Sedimentation(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,input)
      % setup the solver with custom input

      if Sedimentation.checkInput(input)
        error("Simulation not well defined!");
      end

      obj.nmat = length(obj.materials.db)-1;
      obj.grid = gridForSedimentation("XML",input.domain, ...
        "NumMateriais",obj.nmat);

      % Initialize Mesh.
      prepareMesh(obj);

      % Initialize the sedimentation control
      obj.sedimentHistory = SedimentsMap(input.sediment_map,obj.nmat,...
        obj.grid.ncells(1:2));
      obj.sedimentAcc = zeros(prod(obj.grid.ncells(1:2)),obj.nmat);

      % Initialize the BCs
      obj.prepareBC(input.boundary);

      % Initialize the states
      obj.getState().data.(obj.getField()) = zeros(obj.grid.getNumberCells,1);

      % Initialize the Transmissibility.
      prepareSystem(obj);

      obj.fieldId = 1;

      obj.domain.J{obj.fieldId,obj.fieldId} = [];
      obj.domain.rhs{obj.fieldId} = [];

      % Compute the initial system.
      computeHPInitial(obj);

      % Creating the output format

      prepareOutput(obj,input.output);

    end

    function assembleSystem(obj,dt)
      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);
    end

    function advanceState(obj)   % <---- HERE CONTROL TO GROW THE MESH
      % advance the state after reaching convergence

      % Update the sediments accumulated.
      cellGrow = obj.updateSedAccumulated;

      % Update the mesh.


      % hard copy the new state object
      obj.domain.stateOld = copy(obj.domain.state);
    end


    function applyBC(obj,bcId,t)
      if ~BCapplies(obj,bcId)
        return
      end

      [bcDofs,bcVals] = getBC(obj,bcId,t);

      % Base application of a Boundary condition
      bcType = obj.bcs.getType(bcId);

      switch bcType
        case {'Dirichlet','Seepage'}
          applyDirBC(obj,bcId,bcDofs,bcVals);
        case {'Neumann','VolumeForce'}
          applyNeuBC(obj,bcId,bcDofs,bcVals);
        otherwise
          error("Error in %s: Boundary condition type '%s' is not " + ...
            "available in applyBC()",class(obj),bcType)
      end
    end

    function applyDirVal(obj,bcId,t)
      bcVar = obj.bcs.getVariable(bcId);
      if ~strcmp(bcVar,obj.getField())
        return
      end
      [bcDofs,bcVals] = getBC(obj,bcId,t);
      if size(bcVals,2)==2
        % skip BC assigned to external surfaces
        return
      end
      state = getState(obj);
      state.data.pressure(bcDofs) = bcVals;
    end

    % TODO: descomentar
    function updateState(obj,solution)
      ents = obj.grid.getActiveCells;
      state = getState(obj);
      state.data.pressure(ents) = state.data.pressure(ents) + solution;
      % state = getState(obj);
      % state.data.pressure = state.data.pressure + solution;
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      % append state variable to output structure
      sOld = getStateOld(obj);
      sNew = getState(obj);

      p = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);

      outPrint = finalizeState(obj,p,t);
      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma>0
        coords = obj.grid.getCoordCenter(obj.grid.getActiveCells);
        states.potential = p + gamma*coords(:,3);
        states.head = coords(:,3)+p/gamma;
      end
      % mob = (1/obj.materials.getFluid().getDynViscosity());
      % states.flux = computeFlux(obj,p,mob,t);
      states.perm = getPerm(obj,obj.grid.getActiveCells);
      states.pressure = p;
    end

    function writeMatFile(obj,fac,tID)
      pOld = getStateOld(obj,obj.getField());
      pCurr = getState(obj,obj.getField());
      obj.domain.outstate.results(tID).pressure = pCurr*fac+pOld*(1-fac);
    end

    function [dof,vals] = getBC(obj,id,t)
      % getBC - function to find the value and the location for the
      % boundary condition.
      %
      % Observation.:
      %  - The seepage boundary condition apply a hydrostatic pressure
      % in the boundary, and it's assume as a datum the most elavated
      % point in the domain. (For future, have a way to pass this
      % information).
      bc = obj.bcs.db(id);
      [cellId, faceArea, dz] = obj.grid.getBordCell(bc.surface);
      p=obj.domain.state.data.pressure(cellId);
      switch lower(bc.surface)
        case {"x0","xm"}
          axis=1;
          vecN = [1 0 0];
        case {"y0","ym"}
          axis=2;
          vecN = [0 1 0];
        case {"z0","zm"}
          axis=3;
          vecN = [0 0 1];
        otherwise
          dof = [];
          vals = [];
          return
      end

      switch lower(bc.type)
        case 'neumann'
          v=bc.data.*ones(length(cellId),1);
          vals = sum(vecN.*v,2);
        case 'dirichlet'
          gamma = obj.materials.getFluid().getFluidSpecWeight();
          mu = obj.materials.getFluid().getDynViscosity();
          permCell = getPerm(obj,cellId);
          dirJ = 1/mu*(faceArea.*permCell(:,axis));
          potential = p - bc.data + gamma*dz;
          q = dirJ.*potential;
          vals = [dirJ,q];
          % case 'hydrostatic'
          %   gamma = obj.materials.getFluid().getFluidSpecWeight();
          %   zbc = obj.faces.faceCentroid(faceID,3);
          %   href = v(1);
          %   v = gamma*(href-zbc);
          %   v(v<=0)=0.;
          %   mu = obj.materials.getFluid().getDynViscosity();
          %   tr = obj.trans(faceID);
          %   dz = obj.grid.cellCentroid(cellId,3) - obj.faces.faceCentroid(faceID,3);
          %   q = 1/mu*tr.*(p(cellId) - v) + gamma*dz;
          %   vals = [1/mu*tr,q];
      end
      dof = cellId;
    end

    function J = computeMat(obj,dt)
      % recompute elementary matrices only if the model is non-linear
      if obj.grow
        mu = obj.materials.getFluid().getDynViscosity();
        computeStiffMat(obj,1/mu);
        computeCapMat(obj);
      end

      if obj.domain.simparams.isTimeDependent
        J = obj.domain.simparams.theta*obj.H + obj.P/dt;
      else
        J = obj.H;
      end
    end

    function rhs = computeRhs(obj,dt)
      % Compute the residual of the flow problem

      % get pressure state
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());

      lw = 1/obj.materials.getFluid().getDynViscosity();
      ents = obj.grid.mapCellIds(:);
      % ents = obj.dofm.getActiveEntities(obj.fieldId);

      if ~obj.domain.simparams.isTimeDependent
        rhs = obj.H*p(ents);
      else
        theta = obj.domain.simparams.theta;
        rhsStiff = theta*obj.H*p(ents) + (1-theta)*obj.H*pOld(ents);
        rhsCap = (obj.P/dt)*(p(ents) - pOld(ents));
        rhs = rhsStiff + rhsCap;
      end

      %adding gravity rhs contribute
      gamma = obj.materials.getFluid().getFluidSpecWeight();
      if gamma > 0
        rhs = rhs + finalizeRHSGravTerm(obj,lw);
      end

    end

    function computeInitialStiffMat(obj,lw)
      % Compute the Stiffness Matrix for the initial grid
      ncells = obj.grid.numberActiveCells;
      tmpVec = lw.*obj.trans;
      sumDiagTrans = accumarray(obj.facesNeigh(:),repmat(tmpVec,[2,1]),[ncells,1]);
      obj.H = sparse([obj.facesNeigh(:,1); obj.facesNeigh(:,2); (1:ncells)'],...
        [obj.facesNeigh(:,2); obj.facesNeigh(:,1); (1:ncells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], ncells, ncells);
    end

    function computeInitialCapMat(obj)
      ncells = obj.grid.numberActiveCells;
      beta = obj.materials.getFluid().getFluidCompressibility();
      PVal = obj.rockCmCell+beta*obj.poroCell;
      obj.P = PVal.*speye(ncells);

      % subCells = obj.dofm.getFieldCells(obj.fieldId);
      % nSubCells = length(subCells);
      % poroMat = zeros(nSubCells,1);
      % alphaMat = zeros(nSubCells,1);
      % beta = obj.materials.getFluid().getFluidCompressibility();
      % for m = 1:obj.grid.nCellTag
      %   if ~ismember(m,obj.dofm.getTargetRegions([obj.getField(),"displacements"]))
      %     % compute alpha only if there's no coupling in the
      %     % subdomain
      %     alphaMat(m) = obj.materials.getMaterial(m).ConstLaw.getRockCompressibility();
      %   end
      %   poroMat(m) = obj.materials.getMaterial(m).PorousRock.getPorosity();
      % end
      % % (alpha+poro*beta)
      % PVal = alphaMat(obj.grid.cellTag(subCells)) + beta*poroMat(obj.grid.cellTag(subCells));
      % PVal = PVal.*obj.grid.cellVolume(subCells);
      % nDoF = obj.dofm.getNumbDoF(obj.fieldId);
      % [~,~,dof] = unique(subCells);
      % obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
    end

    function gTerm = finalizeRHSGravTerm(obj,lw)
      gTerm = accumarray(obj.facesNeigh(:), ...
        [lw.*obj.rhsGrav; -lw.*obj.rhsGrav],[obj.grid.numberActiveCells,1]);
      gTerm = gTerm(obj.grid.getActiveCells);
    end



    function applyNeuBC(obj,bcId,bcDofs,bcVals)
      if ~BCapplies(obj,bcId)
        return
      end

      % Base application of a Boundary condition
      if ~strcmp(obj.bcs.getVariable(bcId),obj.getField)
        return
      end

      % Base application of Neumann boundary condition to the rhs.
      % bc values are subtracted since we solve du = J\(-rhs)
      bcId = obj.fieldId;
      obj.domain.rhs{bcId}(bcDofs) = obj.domain.rhs{bcId}(bcDofs) - bcVals;
    end

    function applyDirBC(obj,~,bcDofs,bcVals)
      % apply Dirichlet BCs
      % overrides the base method implemented in PhysicsSolver
      % ents: id of constrained faces without any dof mapping applied
      % vals(:,1): Jacobian BC contrib vals(:,2): rhs BC contrib

      assert(size(bcVals,2)==2,'Invalid matrix size for BC values');
      nDoF = obj.grid.numberActiveCells;
      bcDofsJ = nDoF*(bcDofs-1) + bcDofs;
      obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) = ...
        obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) + bcVals(:,1);
      obj.domain.rhs{obj.fieldId}(bcDofs) = ...
        obj.domain.rhs{obj.fieldId}(bcDofs) + bcVals(:,2);
    end

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      pointStr = [];
      % pointStr = repmat(struct('name', 1, 'data', 1), 1, 1);
      % pointStr(1).name = 'flux';
      % pointStr(1).data = state.flux;

      cellStr = repmat(struct('name', 1, 'data', 1), 2, 1);
      cellStr(1).name = 'pressure';
      cellStr(1).data = state.pressure;
      cellStr(2).name = 'permeability';
      cellStr(2).data = state.perm;
      if isfield(state,"potential")
        cellStr(3).name = 'potential';
        cellStr(3).data = state.potential;
        cellStr(4).name = 'piezometric head';
        cellStr(4).data = state.head;
      end
    end

    function out = isLinear(obj)
      out = true;
    end
  end

  methods  (Access = private)

    function prepareMesh(obj)
      lnk = obj.grid;

      obj.mesh.nDim = 3;
      obj.mesh.nCells = lnk.getNumberCells;
      obj.mesh.nNodes = lnk.grid.getNumberPoints;
      obj.mesh.cellVTKType = 12*ones(obj.mesh.nCells,1);
      obj.mesh.cellNumVerts = 8*ones(obj.mesh.nCells,1);
      obj.mesh.cellTag = 8*ones(obj.mesh.nCells,1);

      [obj.mesh.coordinates,obj.mesh.cells] = lnk.grid.getMesh(1);
      obj.mesh.meshType = "Unstructured";
    end

    function prepareBC(obj,data)

      if ~isfield(data,"BC")
        error("No boundary was defined for the simulation");
      end
      nbcs=length(data.BC);
      bc = struct('data',[],'cond','SurfBC',...
        'type',[],'variable',[],'surface',[]);

      for i=1:nbcs
        lnk = data.BC(i);
        bc.type = lnk.type;
        bc.variable = lnk.variable;
        bc.surface = lnk.surface;
        if ~isfield(lnk,'value')
          bc.data=0;
        else
          if isa(lnk.value,"string")
            bc.data=str2num(lnk.value);
          elseif isa(lnk.value,"double")
            bc.data=lnk.value;
          else
            bc.data=0;
          end
        end
        obj.bcs.db(data.BC(i).name) = bc;
      end
    end

    function prepareSystem(obj)
      % Function to create transmissibility for the initial grid
      segmX = diff(obj.grid.grid.X);
      segmY = diff(obj.grid.grid.Y);
      segmZ = diff(obj.grid.grid.Z);

      ncells = prod(obj.grid.ncells);
      cellNeigh = zeros(ncells,6,'uint64');
      cellAreas = zeros(ncells,3);
      for k = 1:obj.grid.ncells(3)
        for i = 1:obj.grid.ncells(1)
          for j = 1:obj.grid.ncells(2)
            cellId = obj.grid.mapCellIds(i, j, k);
            % Neighbor Mapping:
            if i > 1
              cellNeigh(cellId, 3) = obj.grid.mapCellIds(i-1,j,k);
            end
            if i < obj.grid.ncells(1)
              cellNeigh(cellId, 4) = obj.grid.mapCellIds(i+1,j,k);
            end
            if j > 1
              cellNeigh(cellId, 1) = obj.grid.mapCellIds(i,j-1,k);
            end
            if j < obj.grid.ncells(2)
              cellNeigh(cellId, 2) = obj.grid.mapCellIds(i,j+1,k);
            end
            if k > 1
              cellNeigh(cellId, 5) = obj.grid.mapCellIds(i,j,k-1);
            end
            if k < obj.grid.ncells(3)
              cellNeigh(cellId, 6) = obj.grid.mapCellIds(i,j,k+1);
            end

            % Face Area
            cellAreas(cellId, 1) = segmY(j) * segmZ(k);
            cellAreas(cellId, 2) = segmX(i) * segmZ(k);
            cellAreas(cellId, 3) = segmX(i) * segmY(j);
          end
        end
      end

      % Computing the permeability
      permCell = getPerm(obj,obj.grid.getActiveCells);
      ActiveCells = obj.grid.mapCellIds ~= 0;

      % correcting matfrac, cleaning fractions for non-existent cells
      obj.grid.matfrac(~ActiveCells(:),:)=0.;

      % temporay variables
      maxFaces = (obj.grid.ncells(1)+1)*(obj.grid.ncells(2)+1)*(obj.grid.ncells(3)+1);
      tmpNeigh = zeros(maxFaces,2,'uint64');
      isActive = zeros(maxFaces,1,'logical');
      tmpTrans = zeros(maxFaces,1);
      tmpRhs = zeros(maxFaces,1);
      gamma = obj.materials.getFluid().getFluidSpecWeight();

      % Computing transmissibilities and RHS related to the gravity term
      count = 1;
      for cell=1:ncells
        [~,~,zNeiA] = obj.grid.getIJKfromCellID(cell);
        zNeiA=obj.grid.grid.centerZ(zNeiA);
        for face=1:6
          neighbor = cellNeigh(cell, face);
          if (neighbor ~= 0) && (cell < neighbor)
            % Determine axis (1=X, 2=Y, 3=Z)
            axis = mod(ceil(face/2)-1,3)+1;

            % Resistivity: R = 1 / (K * A)
            Ref = 1 / (cellAreas(cell, axis) * permCell(cell, axis));
            Relf = 1 / (cellAreas(neighbor, axis) * permCell(neighbor, axis));

            % Harmonic Average Transmissibility: T = 1 / (Ref + Relf)
            tmpTrans(count) = 1 / (Ref + Relf);

            % Gravity contribution.
            if gamma > 0
              [~,~,zNeiB] = obj.grid.getIJKfromCellID(neighbor);
              zNeiB=obj.grid.grid.centerZ(zNeiB);
              tmpRhs(count) = gamma*tmpTrans(count)*(zNeiA-zNeiB);
            end

            tmpNeigh(count, :) = [cell neighbor];
            isActive(count) = true;
            count = count + 1;
          end
        end
      end

      obj.facesNeigh=tmpNeigh(isActive,:);
      obj.trans=tmpTrans(isActive);
      obj.rhsGrav=tmpRhs(isActive);
    end

    function prepareOutput(obj,data)
      if ~isfield(data,"file")
        error("Output file not pass for the simulation!");
      end

      tmp = OutState(obj.mesh,data.file);

      obj.domain.outstate.modTime = tmp.modTime;
      obj.domain.outstate.timeList = tmp.timeList;
      obj.domain.outstate.results = tmp.results;
      obj.domain.outstate.writeSolution = tmp.writeSolution;
      obj.domain.outstate.timeID = tmp.timeID;
      obj.domain.outstate.writeVtk = tmp.writeVtk;
      obj.domain.outstate.matFileName = tmp.matFileName;
      obj.domain.outstate.vtkFileName = tmp.vtkFileName;

      obj.domain.outstate.VTK = VTKOutput(obj.mesh,obj.domain.outstate.vtkFileName);
    end

    function computeHPInitial(obj)
      ncells = obj.grid.numberActiveCells;
      mu = obj.materials.getFluid().getDynViscosity();

      % Computing the initial H matrix
      tmpVec = (1/mu).*obj.trans;
      sumDiagTrans = accumarray(obj.facesNeigh(:),repmat(tmpVec,[2,1]),[ncells,1]);
      obj.H = sparse([obj.facesNeigh(:,1); obj.facesNeigh(:,2); (1:ncells)'],...
        [obj.facesNeigh(:,2); obj.facesNeigh(:,1); (1:ncells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], ncells, ncells);

      % Computing the initial P matrix
      maxcells = prod(obj.grid.ncells);
      alpha=zeros(maxcells,1);
      poro=zeros(maxcells,1);
      for mats=1:obj.nmat
        tmpMat=obj.materials.getMaterial(mats).ConstLaw.getRockCompressibility();
        alpha=alpha+tmpMat.*obj.grid.matfrac(:,mats);
        tmpMat=obj.materials.getMaterial(mats).PorousRock.getPorosity();
        poro=poro+tmpMat.*obj.grid.matfrac(:,mats);
      end

      % Mapping and Reorder
      mapCellsOn = obj.grid.mapCellIds ~= 0;
      map = obj.grid.mapCellIds(mapCellsOn);
      alpha=alpha(mapCellsOn);
      alpha=alpha(map);
      poro=poro(mapCellsOn);
      poro=poro(map);
      volsCell = obj.grid.computeVols();

      beta = obj.materials.getFluid().getFluidCompressibility();
      PVal = (alpha+beta*poro).*volsCell;
      obj.P = PVal.*speye(ncells);
    end

    function permCells = getPerm(obj,cellId)
      % nelm=length(cellId);
      permCells=0;
      for mat=1:obj.nmat
        tmpMat=obj.materials.getMaterial(mat).PorousRock.getPermVoigt();
        permCells = permCells+obj.grid.matfrac(cellId,mat)*tmpMat(1:3);
      end
    end

    function cellGrow = updateSedAccumulated(obj)
      % UPDATESEDACCUMULATED Updates sediment buffer and triggers grid growth.
      %
      %   Calculates new deposition, updates the accumulation buffer, and
      %   identifies columns reaching the height threshold for new cells.
      %
      %   Returns:
      %       cellGrow - Logical mask of columns requiring a new cell layer.

      % 1. Calculate time step
      t0 = obj.domain.stateOld.t;
      dt = obj.domain.state.t-t0;

      % 2. Update accumulation buffer
      addSed = obj.sedimentHistory.getSedimentationMap(t0,dt);
      obj.sedimentAcc = obj.sedimentAcc + addSed;

      % 3. Check for growth trigger (Height >= Threshold)
      colSed = sum(obj.sedimentAcc,2);
      cellGrow = colSed >= obj.heighControl;

      % 4. Handle overflow for grown cells
      if any(cellGrow)
        dl = obj.heighControl-colSed;
        sed = sum(addSed,2);
        obj.sedimentAcc(cellGrow,:)=(1/dt)*(dl(cellGrow)./sed(cellGrow)).*addSed(cellGrow,:);
      end
    end



  end



  methods (Static)

    function out = getField()
      out = "pressure";
    end

    function out = isSymmetric()
      out = true;
    end

    function flag = checkInput(input)
      flag = false;
      if ~(isfield(input,'domain') || isfield(input,'Domain'))
        flag = true;
        disp("The initial domain for the simulation is not defined!");
      end

      if ~(isfield(input,'boundary') || isfield(input,'Boundary'))
        flag = true;
        disp("The boundary condition for the simulation is not defined!");
      end

      if ~(isfield(input,'sediment_map') || isfield(input,'Sediment_map'))
        flag = true;
        disp("Map of sedimentation for your simulation not defined!");
      end

      if ~(isfield(input,'output') || isfield(input,'Output'))
        flag = true;
        disp("The output specifications is not defined!");
      end
    end

  end

end


