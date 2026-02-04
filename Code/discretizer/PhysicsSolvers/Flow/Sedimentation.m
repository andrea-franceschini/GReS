classdef Sedimentation < PhysicsSolver
  % SEDIMENTATION
  % ------------------------------------------------------------------
  % Transient single-phase Darcy flow solver with sediment-driven
  % vertical mesh growth.
  %
  % Governing physics:
  %   - Darcy flow with gravity
  %   - Rock and fluid compressibility
  %
  % Discretization:
  %   - Finite volume
  %   - Theta-method time integration
  %
  % Mesh:
  %   - Structured column-based grid
  %   - Dynamic growth in Z direction due to sediment accumulation
  %
  % Growth criterion:
  %   - A new cell is added when accumulated sediment height
  %     exceeds heightControl.
  %
  % Notes:
  %   - Grid topology, transmissibility, and state vectors
  %     may change during the simulation.
  % ------------------------------------------------------------------

  properties
    grid gridForSedimentation           % Grid with column growth support

    sedimentHistory = SedimentsMap()    % Time-dependent sediment control
    sedimentAcc (:,:)                   % Accumulated sediment per column
    heightControl double = 0.1          % Height threshold for new cell

    facesNeigh (:,2)                    % Neighboring cell pairs per face
    facesNeighDir (:,1)                 % Face normal direction (1=x,2=y,3=z)
    halfTrans (:,3)                     % Half transmissibility per cell

    H (:,:)
    P (:,:)
    void0 (:,1) double                  % Initial porosity for each cell
  end

  properties (Access = protected)
    fieldId
  end

  properties (Access =private)
    nmat uint16
    sedFlag logical = true
  end

  methods (Access = public)
    function obj = Sedimentation(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,input)
      % REGISTERSOLVER Initializes solver data structures.

      % Validate input
      if Sedimentation.checkInput(input)
        gresLog().error("Simulation not well defined!");
      end

      % Build grid & mesh
      obj.nmat = length(obj.materials.db)-1;
      obj.grid = gridForSedimentation("XML",input.Domain, ...
        "NumMateriais",obj.nmat);
      if isfield(input.Domain,"NewCellHeightControl")
        obj.heightControl = input.Domain.NewCellHeightControl;
      end
      obj.prepareMesh;

      % Initialize sediment control
      obj.sedimentHistory = SedimentsMap(input.SedimentMap,obj.nmat,...
        obj.grid.ncells(1:2));
      obj.sedimentAcc = zeros(prod(obj.grid.ncells(1:2)),obj.nmat);

      % Setup BCs
      obj.prepareBC(input.Boundary);

      % Initialize the States
      obj.prepareStates;

      % Allocate system matrices
      obj.prepareSystem;
      obj.fieldId = 1;
      obj.domain.J{obj.fieldId,obj.fieldId} = [];
      obj.domain.rhs{obj.fieldId} = [];

      % Increase number of variables in the dof manager
      obj.dofm.registerVariable(obj.getField(),entityField.cell,1);

      % Prepare output
      obj.prepareOutput(input.Output);
    end

    function assembleSystem(obj,dt)
      % Only one time for each time step, update the sedimentation rate.
      if obj.sedFlag
        obj.updateSedRate(dt);
        obj.sedFlag = false;
      end

      % Assembling the system.
      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);
    end

    function advanceState(obj)
      % ADVANCESTATE Finalizes time step and updates grid topology.

      % Update sediment accumulation
      [cellGrow, cellSed] = obj.updateSedAccumulated;

      % Grow mesh if height threshold is reached
      flagGrow = any(cellGrow);
      if flagGrow
        meshUpdate(obj,cellGrow,cellSed);
      end

      % Update state for next step
      obj.domain.stateOld = copy(obj.domain.state);

      % Update the dof used for the parallel solver
      if flagGrow
        obj.getState().data.maxDofUnchanged = obj.grid.getMaxDofUnchanged;
      end

      % Update the necessity to compute the sedimentation rate.
      obj.sedFlag = true;
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
          gresLog().error("Error in %s: Boundary condition type '%s' is not " + ...
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

    function updateState(obj,solution)
      state = getState(obj);
      state.data.pressure = state.data.pressure + solution;
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      % append state variable to output structure
      % sOld = getStateOld(obj);
      % sNew = getState(obj);
      % p = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);

      sOld = obj.getStateOld(obj.getField);
      sNew = obj.getState(obj.getField);
      p = sNew*fac+sOld*(1-fac);

      outPrint = finalizeState(obj,p,t);
      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function states = finalizeState(obj,p,t)
      % Compute the posprocessing variables for the module.
      gamma = obj.materials.getFluid().getSpecificWeight();
      if gamma>0
        coords = obj.grid.getCoordCenter(obj.grid.getActiveDofs);
        states.potential = p + gamma*coords(:,3);
        states.head = coords(:,3)+p/gamma;
      end
      % mob = (1/obj.materials.getFluid().getDynViscosity());
      % states.flux = computeFlux(obj,p,mob,t);

      states.poro = obj.grid.accumulateProp(obj.getProp('porosity'));
      states.perm = obj.grid.accumulateProp(obj.getProp('permeability'));
      states.comp = obj.getState().data.cellComp;
      states.stress = obj.getState().data.stress;
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
      p = getState(obj,"pressure");
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
          gamma = obj.materials.getFluid().getSpecificWeight();
          mu = obj.materials.getFluid().getDynViscosity();
          permCell = obj.grid.accumulateProp(obj.getProp('permeability'),cellId);
          dirJ = 1/mu*(faceArea.*permCell(:,axis));
          % potential = p(cellId) - bc.data - gamma*dz;
          potential = p(cellId) - bc.data;
          q = dirJ.*potential;
          vals = [dirJ,q];
      end
      dof = cellId;
    end

    function J = computeMat(obj,dt)
      % Recompute elementary matrices
      obj.computeStiffMat;
      obj.computeCapMat;
      J = obj.domain.simparams.theta*obj.H + obj.P/dt;
    end

    function rhs = computeRhs(obj,dt)
      % Compute the residual of the flow problem

      % get pressure state
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());

      % lw = 1/obj.materials.getFluid().getDynViscosity();
      theta = obj.domain.simparams.theta;
      rhsStiff = theta*obj.H*p + (1-theta)*obj.H*pOld;
      rhsCap = (obj.P/dt)*(p - pOld);
      rhs = rhsStiff + rhsCap;

      % adding sediment contribution
      [idI, idJ, idK] = obj.grid.getIJKTop;
      dofs = obj.grid.getDofsFromIJK([idI,idJ,idK]);
      rhs(dofs) = rhs(dofs) + ...
        obj.computeOedometricCompressibility(dofs).*obj.getState('tstressvar');

      % % % %adding gravity rhs contribute
      % % % gamma = obj.materials.getFluid().getSpecificWeight();
      % % % if gamma > 0
      % % %   rhs = rhs + computeGravTerm(obj,lw);
      % % % end
    end

    function computeStiffMat(obj)
      % Compute the Stiffness Matrix
      lw = 1/obj.materials.getFluid().getDynViscosity();
      ncells = obj.grid.ndofs;
      tmpVec = lw.*obj.computeTrans;
      sumDiagTrans = accumarray(obj.facesNeigh(:),repmat(tmpVec,[2,1]),[ncells,1]);
      obj.H = sparse([obj.facesNeigh(:,1); obj.facesNeigh(:,2); (1:ncells)'],...
        [obj.facesNeigh(:,2); obj.facesNeigh(:,1); (1:ncells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], ncells, ncells);
    end

    function computeCapMat(obj)
      % COMPUTECAPMAT Builds the storage matrix (P).
      %
      % Includes:
      %   - Rock compressibility
      %   - Fluid compressibility

      % Selecting the active dofs
      % dofs = obj.grid.getActiveDofs;
      dofs = 1:obj.grid.ndofs;

      % Computing the volumes
      volsCell = obj.grid.computeVols(dofs,obj.getState().data.cellComp);

      % Computing the storage coefficient
      poro = obj.getState().data.voidrate;
      poro = poro./(1+poro);
      beta = obj.materials.getFluid().getFluidCompressibility();
      oedoComp = obj.computeOedometricCompressibility(dofs);
      PVal = (oedoComp+beta*poro).*volsCell;

      % Creating the P matrix.
      obj.P = PVal.*speye(obj.grid.ndofs);
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
      nDoF = obj.grid.ndofs;
      bcDofsJ = nDoF*(bcDofs-1) + bcDofs;
      obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) = ...
        obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) + bcVals(:,1);
      obj.domain.rhs{obj.fieldId}(bcDofs) = ...
        obj.domain.rhs{obj.fieldId}(bcDofs) + bcVals(:,2);
    end

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      pointStr = [];

      cellStr = repmat(struct('name', 1, 'data', 1), 2, 1);
      cellStr(1).name = 'pressure';
      cellStr(1).data = state.pressure;
      cellStr(2).name = 'permeability';
      cellStr(2).data = state.perm;
      cellStr(3).name = 'porosity';
      cellStr(3).data = state.poro;
      cellStr(4).name = 'Compaction(cell)';
      cellStr(4).data = state.comp;
      cellStr(5).name = 'stress(vertical)';
      cellStr(5).data = state.stress;
      if isfield(state,"potential")
        cellStr(6).name = 'potential';
        cellStr(6).data = state.potential;
        cellStr(7).name = 'piezometric head';
        cellStr(7).data = state.head;
      end
    end

    function out = isLinear(obj)
      out = true;
    end
  end

  methods (Access = private)

    function prepareMesh(obj)
      obj.mesh.nDim = 3;
      obj.mesh.nCells = obj.grid.getNumberCells;
      obj.mesh.nNodes = obj.grid.getNumberPoints;
      obj.mesh.cellVTKType = 12*ones(obj.mesh.nCells,1);
      obj.mesh.cellNumVerts = 8*ones(obj.mesh.nCells,1);
      obj.mesh.cellTag = 8*ones(obj.mesh.nCells,1);

      [obj.mesh.coordinates,obj.mesh.cells] = obj.grid.getMesh;
      obj.mesh.meshType = "Unstructured";
    end

    function prepareBC(obj,data)

      if ~isfield(data,"BC")
        gresLog().error("No boundary was defined for the simulation");
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

    function prepareStates(obj)
      ndofs = obj.grid.getNumberCells;
      obj.getState().data.maxDofUnchanged = 0;
      obj.getState().data.(obj.getField()) = zeros(ndofs,1);
      obj.getState().data.cellComp = zeros(ndofs,1);
      obj.getState().data.stress = zeros(ndofs,1);
      obj.getState().data.sedmrate = zeros(ndofs,1);
      obj.getState().data.tstressvar = [];
      obj.getState().data.voidrate = zeros(ndofs,1);
    end

    function prepareSystem(obj)
      % Compute the Neighborhood.
      dofs = obj.grid.getDofsFromIJK;
      ActiveCells = dofs ~= 0;

      actDofs = dofs(ActiveCells);
      actIJK = obj.grid.getIJKfromCellID(actDofs);
      cellNeigh = obj.grid.getNeigh(actIJK);
      celldims = obj.grid.getDims(actIJK);

      newcells=sum(ActiveCells);
      cellsConc = zeros([6*newcells,2]);
      cellsConcDir = ones([6*newcells,1]);
      cellsConcActive = false([6*newcells,1]);
      for ref=1:6
        % Add the connection between the reference cell and the neighborhoods
        tmp = cellNeigh(:,ref);

        % Evaluate the neighbor is smaller than the dof.
        flag = tmp > actDofs;
        nfaces = sum(flag);

        cellsConc(ref*newcells+1:ref*newcells+nfaces,:)=[actDofs(flag),tmp(flag)];
        cellsConcDir(ref*newcells+1:ref*newcells+nfaces)=ceil(ref/2);
        cellsConcActive(ref*newcells+1:ref*newcells+nfaces)=true;
      end
      obj.facesNeigh = cellsConc(cellsConcActive,:);
      obj.facesNeighDir = cellsConcDir(cellsConcActive);

      % Computing the permeability
      permCell = obj.grid.accumulateProp(obj.getProp('permeability'));

      % correcting matfrac, cleaning fractions for non-existent cells
      % obj.grid.matfrac(~ActiveCells(:),:)=0.;

      % Computing half transmissibilities
      obj.halfTrans = zeros(obj.grid.ndofs,3);
      obj.halfTrans(:,1)=1./(celldims(actDofs,2).*permCell(actDofs,1));
      obj.halfTrans(:,2)=1./(celldims(actDofs,1).*permCell(actDofs,2));
      obj.halfTrans(:,3)=1./(celldims(actDofs,1).*celldims(actDofs,2).*permCell(actDofs,3));

      % Computing the initial void rate
      obj.void0 = obj.grid.accumulateProp(obj.getProp('voidrate'));
      obj.getState().data.voidrate = obj.void0./(1+obj.void0);
    end

    function prepareOutput(obj,data)
      if ~isfield(data,"file")
        gresLog().error("Output file not pass for the simulation!");
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

    function updateSedRate(obj,dt)
      % Finding the dof for the top of each column in the grid
      [idI, idJ, idK] = obj.grid.getIJKTop;
      dofs = obj.grid.getDofsFromIJK([idI,idJ,idK]);

      % Finding the sedimentation rate
      t0 = obj.getState().t;
      sedRate = obj.sedimentHistory.getSedimentationMap(t0,dt);
      obj.getState().data.sedmrate = sedRate;

      % Finding the total stress variation
      poro = 1 - obj.void0(dofs)./(1+obj.void0(dofs));
      spwg = obj.grid.accumulateProp(obj.getProp('specificweight'),dofs) ...
        - obj.materials.getFluid().getSpecificWeight();
      obj.getState().data.tstressvar = poro.*spwg.*sedRate;
    end

    function [cellGrow, cellSed] = updateSedAccumulated(obj)
      % UPDATESEDACCUMULATED Updates sediment buffer and triggers grid growth.
      %
      %   Calculates new deposition, updates the accumulation buffer, and
      %   identifies columns reaching the height threshold for new cells.
      %
      % Outputs:
      %   cellGrow - Logical array indicating columns that grow
      %   cellSed  - Sediment assigned to newly created cells

      % 0. Initializate cellSed
      cellSed = zeros(size(obj.sedimentAcc));
      
      % 2. Update accumulation buffer
      addSed = obj.getState().data.sedmrate;
      obj.sedimentAcc = obj.sedimentAcc + addSed;

      % 3. Check for growth trigger (Height >= Threshold)
      colSed = sum(obj.sedimentAcc,2);
      cellGrow = colSed >= obj.heightControl;

      % 4. Handle overflow for grown cells
      if any(cellGrow)
        dl = colSed-obj.heightControl;
        sed = (dl./sum(addSed,2)).*addSed;
        cellSed(cellGrow,:) = obj.sedimentAcc(cellGrow,:)-sed(cellGrow,:);
        obj.sedimentAcc(cellGrow,:) = sed(cellGrow,:);
      end
    end

    function meshUpdate(obj,map,sed)
      % MESHUPDATE Updates grid topology due to sediment growth.
      %
      % Actions:
      %   - Add new cells
      %   - Update face connectivity
      %   - Update transmissibility
      %   - Extend state vectors
      %   - Update VTK mesh

      lnk=obj.grid;
      newcells = sum(map);

      colheightBG = lnk.columnsHeight;
      matfrac = sed/obj.heightControl;
      colNotTop = colheightBG~=lnk.ncells(3);
      colNotTop = colNotTop(map);

      % Update the grid and locate the new dofs
      newlayer = lnk.grow(map,matfrac,obj.heightControl);
      [idI,idJ,idK] = lnk.getIJKTop;
      actIJK = [idI(map),idJ(map),idK(map)];
      dofs = lnk.getDofsFromIJK(actIJK);

      % Update initial void rate
      voidI = obj.grid.accumulateProp(obj.getProp('voidrate'),dofs);
      obj.void0(end+1:end+newcells) = voidI;

      % Update the face connectivity
      permCell = obj.grid.accumulateProp(obj.getProp('permeability'),dofs);
      celldims = lnk.getDims(actIJK);

      % Creating the connectives between new cells
      cellNeigh = lnk.getNeigh(actIJK);

      cellsConc = zeros([5*newcells,2]);
      cellsConcDir = ones([5*newcells,1]);
      cellsConcActive = false([5*newcells,1]);
      cellsConc(1:newcells,:)=[cellNeigh(:,5) dofs];
      cellsConcDir(1:newcells)=3;
      cellsConcActive(1:newcells)=true;
      for ref=1:4
        % Add the connection between the reference cell and the neighborhoods
        tmp = cellNeigh(:,ref);

        % Evaluate the neighbor if the column is not at the top.
        flag0 = and(and(tmp<dofs,tmp~=0),colNotTop);

        % Evaluate the neighbor if the column is at the top.
        % flag1 = and(tmp > dofs,~colNotTop);
        flag1 = tmp > dofs;

        % Combine the two case.
        flag = or(flag0,flag1);
        nfaces = sum(flag);

        cellsConc(ref*newcells+1:ref*newcells+nfaces,:)=[dofs(flag),tmp(flag)];
        cellsConcDir(ref*newcells+1:ref*newcells+nfaces)=ceil(ref/2);
        cellsConcActive(ref*newcells+1:ref*newcells+nfaces)=true;
      end
      nfaces = sum(cellsConcActive);
      obj.facesNeigh(end+1:end+nfaces,:) = cellsConc(cellsConcActive,:);
      obj.facesNeighDir(end+1:end+nfaces,:) = cellsConcDir(cellsConcActive);

      % Computing the new half transmissibilities
      appAtEnd = size(obj.halfTrans,1);
      obj.halfTrans(appAtEnd+1:appAtEnd+newcells,1) = ...
        1./(celldims(:,2).*permCell(:,1));
      obj.halfTrans(appAtEnd+1:appAtEnd+newcells,2) = ...
        1./(celldims(:,1).*permCell(:,2));
      obj.halfTrans(appAtEnd+1:appAtEnd+newcells,3) = ...
        1./(celldims(:,1).*celldims(:,2).*permCell(:,3));

      % Update the states
      % idK = colheightBG;
      % actIJK = [idI(map),idJ(map),idK(map)];
      % dofs = lnk.getDofsFromIJK(actIJK);
      % obj.domain.state.data.pressure(end+1:end+newcells) = ...
      %   obj.domain.state.data.pressure(dofs);

      obj.domain.state.data.pressure(end+1:end+newcells) = 0.;

      % poro = zeros(newcells,1);
      % for mats=1:obj.nmat
      %   pos = sub2ind(size(obj.grid.matfrac), dofs,repelem(mats,newcells,1));
      %   tmpMat=obj.materials.getMaterial(mats).PorousRock.getPorosity();
      %   poro=poro+tmpMat.*obj.grid.matfrac(pos);
      % end

      % Update the cell displacement.
      obj.domain.state.data.cellComp(end+1:end+newcells) = 0.;
      obj.domain.state.data.stress(end+1:end+newcells) = 0.;
      obj.domain.state.data.voidrate(end+1:end+newcells) = voidI;

      % Update the mesh output
      obj.updateMeshOutput(map,newlayer);
    end

    function updateMeshOutput(obj,map,newlayer)
      newcells = sum(map);
      if newlayer
        [XX, YY, ZZ] = ndgrid(obj.grid.coordX, obj.grid.coordY, obj.grid.coordZ(end));
        coord = [XX(:), YY(:), ZZ(:)];
        npoints = size(coord,1);
        obj.mesh.nNodes = obj.mesh.nNodes+npoints;
        obj.mesh.coordinates(end+1:obj.mesh.nNodes,:) = coord;
      end

      [idI,idJ,idK] = obj.grid.getIJKTop;

      obj.mesh.cells(end+1:end+newcells,:) = ...
        obj.grid.getConectByIJK(idI(map),idJ(map),idK(map));
      obj.mesh.nCells=obj.mesh.nCells+newcells;
      obj.mesh.cellTag(end+1:end+newcells) = 8;
      obj.mesh.cellNumVerts(end+1:end+newcells) = 8;
      obj.mesh.cellVTKType(end+1:end+newcells) = 12;
    end

    

    


    

    function transm = computeTrans(obj)
      % COMPUTETRANS Computes face transmissibilities.
      %
      % Notes:
      %   - X and Y half transmissibilities are stored without height.

      lnk=obj.grid;

      % Computing the permeability
      actDof = lnk.dof ~=0;
      [idI,idJ,idK] = lnk.getIJK;

      % Finalize the half-transmissibility
      dims = lnk.getDims([idI(actDof),idJ(actDof),idK(actDof)]);
      tmpT = obj.halfTrans;
      dl = obj.getState().data.cellComp;
      for i=1:2
        tmpT(:,i) = tmpT(:,i)./(dims(:,3)-dl);
        % tmpT(:,i) = tmpT(:,i)./dims(:,3);   % <--- without cell height update
      end

      idx = sub2ind(size(tmpT), obj.facesNeigh(:,1), obj.facesNeighDir);
      Tii = tmpT(idx);
      idx = sub2ind(size(tmpT), obj.facesNeigh(:,2), obj.facesNeighDir);
      Tik = tmpT(idx);
      transm = 1./(Tii+Tik);
    end

    function gTerm = computeGravTerm(obj,lw)
      % COMPUTEGRAVTERM Computes gravity contribution to RHS.

      % Computing the permeability
      dofs = obj.grid.getActiveDofs;
      ijk = obj.grid.getIJKfromCellID(dofs);
      zcells = obj.grid.coordZ(1:end-1)+diff(obj.grid.coordZ)/2.;
      % zcells = zcells(ijk(:,3));
      zcells = zcells(ijk(:,3))-obj.domain.state.data.cellComp; % update the z height
      zneiA = zcells(obj.facesNeigh(:,1));
      zneiB = zcells(obj.facesNeigh(:,2));

      gamma = obj.materials.getFluid().getSpecificWeight();
      tmpVec = gamma*lw.*obj.computeTrans.*(zneiA-zneiB);
      gTerm = accumarray(obj.facesNeigh(:), ...
        [tmpVec; -tmpVec],[obj.grid.ndofs,1]);
    end

    function oedoComp = computeOedometricCompressibility(obj,dofs)
      if ~exist("dofs","var")
        dofs = obj.grid.getActiveDofs;
      end
      oedoComp=zeros(length(dofs),1);
    end






    function out = getProp(obj,type)
      switch type
        case 'porosity'
          out = zeros(obj.nmat,1);
          for mat=1:obj.nmat
            out(mat)=obj.materials.getMaterial(mat).PorousRock.getPorosity();
          end
        case 'specificweight'
          out = zeros(obj.nmat,1);
          for mat=1:obj.nmat
            out(mat)=obj.materials.getMaterial(mat).PorousRock.getSpecificWeight();
          end
        case 'permeability'
          out = zeros(obj.nmat,3);
          for mat=1:obj.nmat
            out(mat,:)=diag(obj.materials.getMaterial(mat).PorousRock.getPermMatrix());
          end
        case 'voidrate'
          out = zeros(obj.nmat,1);
          for mat=1:obj.nmat
            out(mat,:)=diag(obj.materials.getMaterial(mat).SedMaterial.getVoidRate());
          end
        case 'compressIdx'
          out = zeros(obj.nmat,1);
          for mat=1:obj.nmat
            out(mat,:)=diag(obj.materials.getMaterial(mat).SedMaterial.getCompressibilityIdx());
          end
        case 'recompressIdx'
          out = zeros(obj.nmat,1);
          for mat=1:obj.nmat
            out(mat,:)=diag(obj.materials.getMaterial(mat).SedMaterial.getReCompressibilityIdx());
          end
        case 'preConStress'
          out = zeros(obj.nmat,1);
          for mat=1:obj.nmat
            out(mat,:)=diag(obj.materials.getMaterial(mat).SedMaterial.getPreConsolidadeStress());
          end
        otherwise
          out = [];
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
      if ~isfield(input,'Domain')
        flag = true;
        disp("The initial domain for the simulation is not defined!");
      end

      if ~isfield(input,'Boundary')
        flag = true;
        disp("The boundary condition for the simulation is not defined!");
      end

      if ~isfield(input,'SedimentMap')
        flag = true;
        disp("Map of sedimentation for your simulation not defined!");
      end

      if ~isfield(input,'Output')
        flag = true;
        disp("The output specifications is not defined!");
      end
    end


  end

end