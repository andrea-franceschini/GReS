classdef Sedimentation < PhysicsSolver
  % SEDIMENTATION
  % ------------------------------------------------------------------
  % Transient single-phase Darcy flow solver with sediment-driven
  % vertical mesh growth.
  %
  % Governing physics:
  %   - Darcy flow
  %   - Rock and fluid compressibility
  %
  % Discretization:
  %   - Finite volume
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
    
    poro0 (:,1) double                  % Initial porosity for each cell
    matfrac (:,:)                       % Material fractions per cell
  end

  properties (Access = protected)
    fieldId
  end

  properties (Access = private)
    nmat uint16
    sedFlag logical = true
    nonElasticFlag logical = true
  end

  methods (Access = public)
    function obj = Sedimentation(domain)
      obj@PhysicsSolver(domain);
    end

    function registerSolver(obj,input)
      % REGISTERSOLVER Initializes solver data structures.

      % Validate input
      obj.checkInput(input);

      % Build grid
      obj.grid = gridForSedimentation("XML",input.Domain);

      % Build the material fractions for each cell
      obj.prepareMaterialFractions(input.Domain.Initial);

      % Set the default height for news cells.
      if isfield(input.Domain,"NewCellHeightControl")
        obj.heightControl = input.Domain.NewCellHeightControl;
      end

      % Set mesh to be print.
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
    end

    function assembleSystem(obj,dt)
      % Update the cells half transmissibility
      obj.computeHalfTrans();

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
      if obj.nonElasticFlag
        map = obj.domain.state.data.stress > obj.domain.state.data.prestress;
        obj.domain.state.data.prestress(map) = obj.domain.state.data.stress(map);
      end

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

    % function applyDirVal(obj,bcId,t)
    %   bcVar = obj.bcs.getVariable(bcId);
    %   if ~strcmp(bcVar,obj.getField())
    %     return
    %   end
    %   [bcDofs,bcVals] = getBC(obj,bcId,t);
    %   if size(bcVals,2)==2
    %     % skip BC assigned to external surfaces
    %     return
    %   end
    %   state = getState(obj);
    %   state.data.pressure(bcDofs) = bcVals;
    % end
    function applyDirVal(obj,bcId,t)
      return
    end

    function updateState(obj,solution)
      % Update overpressure
      state = getState(obj);
      stateOld = getStateOld(obj);
      state.data.pressure = state.data.pressure + solution;
      dp = state.data.pressure - stateOld.data.pressure;

      % Update the stress state
      dt = state.t-obj.domain.stateOld.t;
      sOld = stateOld.data.stress;

      topCells = obj.grid.getBordCell("ZM");
      gamma_s = obj.getCellsProp('specificWeight',topCells);
      gamma_w = obj.materials.getFluid().getSpecificWeight();
      iniPoro = obj.poro0(topCells);
      dofList = obj.grid.getMapFormDofs();
      omega = obj.getState().data.sedmrate;
      K = (gamma_s - gamma_w) .* (1 - iniPoro) .* omega;
      sNew = sOld + dp - dt*K(dofList);

      % map = reshape(state.data.tstressvar,obj.grid.ncells(1:2));
      % sNew = sOld-dt*obj.grid.distMapOverDofs(map)+dp;
      % sNew = sOld + dp - dt*(1-obj.poro0).*obj.getCellsProp('deltaWeight') ...
      %   .* obj.getState().data.sedmrate(obj.grid.getMapFormDofs());
      state.data.stress = sNew;


      
      % Update the deformation      
      if obj.nonElasticFlag
        obj.updateStateNonElastic(sNew,sOld);
      else
        obj.updateStateElastic(sNew,sOld);
      end
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      % append state variable to output structure
      sOld = obj.getStateOld();
      sNew = obj.getState();
      p = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);

      outPrint.overpres = p;
      % outPrint.overpres = p/gamma;
      outPrint.perm = obj.getCellsProp('conductivity');
      comp = sNew.data.cellDefm*fac+sOld.data.cellDefm*(1-fac);
      outPrint.comp = obj.grid.getCompaction(comp);
      outPrint.stress = sNew.data.stress*fac+sOld.data.stress*(1-fac);

      gamma = obj.materials.getFluid().getSpecificWeight();
      coords = obj.grid.getCoordCenter(obj.grid.getActiveDofs);
      outPrint.pressure = p + gamma*coords(:,3);
      outPrint.head = coords(:,3)+p/gamma;

      if obj.nonElasticFlag
        voidR = sNew.data.voidrate*fac+sOld.data.voidrate*(1-fac);
        outPrint.poro = voidR./(1+voidR);
      else
        outPrint.poro = obj.poro0;
      end

      % mob = (1/obj.materials.getFluid().getDynViscosity());
      % states.flux = computeFlux(obj,p,mob,t);

      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function writeMatFile(obj,fac,tID)
      sOld = obj.getStateOld();
      sNew = obj.getState();

      press = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);
      obj.domain.outstate.matFile(tID).pressure = press;

      comp = sNew.data.cellDefm*fac+sOld.data.cellDefm*(1-fac);
      obj.domain.outstate.matFile(tID).compaction = obj.grid.getCompaction(comp);
      % obj.domain.outstate.matFile(tID).
      % 
      % 
      % outPrint.comp = 
      % outPrint.stress = sNew.data.stress*fac+sOld.data.stress*(1-fac);
      % 
      % gamma = obj.materials.getFluid().getSpecificWeight();
      % coords = obj.grid.getCoordCenter(obj.grid.getActiveDofs);
      % outPrint.pressure = p + gamma*coords(:,3);
      % outPrint.head = coords(:,3)+p/gamma;
      % 
      % if obj.nonElasticFlag
      %   voidR = sNew.data.voidrate*fac+sOld.data.voidrate*(1-fac);
      %   outPrint.poro = voidR./(1+voidR);
      % else
      %   outPrint.poro = obj.poro0;
      % end

      % mob = (1/obj.materials.getFluid().getDynViscosity());
      % states.flux = computeFlux(obj,p,mob,t);

      % [cellData,pointData] = buildPrintStruct(obj,outPrint);

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
      cellId = obj.grid.getBordCell(bc.surface);
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
          mu = obj.materials.getFluid().getSpecificWeight();
          dirJ = 1/mu*obj.halfTrans(cellId,axis);
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

      theta = obj.domain.simparams.theta;
      rhsStiff = theta*obj.H*p + (1-theta)*obj.H*pOld;
      rhsCap = (obj.P/dt)*(p - pOld);


      topCells = obj.grid.getBordCell("ZM");
      gamma_s = obj.getCellsProp('specificWeight',topCells);
      gamma_w = obj.materials.getFluid().getSpecificWeight();
      iniPoro = obj.poro0(topCells);
      dofList = obj.grid.getMapFormDofs();
      omega = obj.getState().data.sedmrate;
      K = (gamma_s - gamma_w) .* (1 - iniPoro) .* omega;
      cb = obj.computeOedometricCompressibility();

      [dx,dy,dz] = obj.grid.getCellsDims();
      volCell = dx.*dy.*(dz + obj.getState('cellDefm'));      
      rhsSedm = volCell .* cb .* K(dofList);


      % rhsSedm = volCell ...
      % rhsSedm = (1/dt)*volCell ...  
      %   .* obj.computeOedometricCompressibility() .* (1-obj.poro0) ...
      %   .* obj.getCellsProp('deltaWeight') .* obj.getState().data.sedmrate(obj.grid.getMapFormDofs());

      rhs = rhsStiff + rhsCap - rhsSedm;

      % adding sediment contribution
      % map = reshape(obj.getState('tstressvar'),obj.grid.ncells(1:2));
      % rhs = rhs - volCell.*obj.computeOedometricCompressibility().*obj.grid.distMapOverDofs(map);
    end

    function computeStiffMat(obj)
      % Compute the Stiffness Matrix
      ncells = obj.grid.ndofs;
      lw = 1/obj.materials.getFluid().getSpecificWeight();

      idi = sub2ind(size(obj.halfTrans), obj.facesNeigh(:,1), obj.facesNeighDir);
      idk = sub2ind(size(obj.halfTrans), obj.facesNeigh(:,2), obj.facesNeighDir);
      Tii = obj.halfTrans(idi);      
      Tik = obj.halfTrans(idk);
      trans = 1./(1./Tii+1./Tik);
      tmpVec = lw.*trans;
      
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
      [dx,dy,dz] = obj.grid.getCellsDims();
      volCell = dx.*dy.*(dz + obj.getState('cellDefm'));

      % Computing the storage coefficient
      if obj.nonElasticFlag
        poro = obj.getState().data.voidrate;
        poro = poro./(1+poro);
      else
        poro = obj.poro0;
      end
      beta = obj.materials.getFluid().getFluidCompressibility();
      oedoComp = obj.computeOedometricCompressibility(dofs);
      PVal = (oedoComp+beta*poro).*volCell;

      % CreatcellDefming the P matrix.
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
      cellStr(1).name = 'overpressure';
      cellStr(1).data = state.overpres;
      cellStr(2).name = 'pressure';
      cellStr(2).data = state.pressure;
      cellStr(3).name = 'piezometric head';
      cellStr(3).data = state.head;

      cellStr(4).name = 'conductivity';
      cellStr(4).data = state.perm;
      cellStr(5).name = 'porosity';
      cellStr(5).data = state.poro;
      cellStr(6).name = 'compaction(cell)';
      cellStr(6).data = -state.comp;
      cellStr(7).name = 'stress(vertical)';
      cellStr(7).data = -state.stress;
    end

    function out = isLinear(obj)
      out = true;
    end
  end

  methods (Access = private)

    function checkInput(obj,input)
      flag = false;

      % Check the material is well defined.
      obj.nmat = length(obj.materials.db)-1;
      if isfield(obj.materials.db(1),'Elastic')
        tmp1=true;
      else
        tmp1=false;
      end      
      for i=1:obj.nmat
        tmp2=false;
        if isfield(obj.materials.db(i),'Elastic')
          tmp2=true;          
        end
        if tmp2 ~= tmp1
          flag = true;
          disp("The materials for the simulation is not well defined!");
        end
        tmp1=tmp2;
      end
      obj.nonElasticFlag = ~tmp1;

      % Check other input      
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

      if flag
        gresLog().error("Simulation is not well defined!");
      end
    end

    function prepareMaterialFractions(obj,data)
      % Build the material fractions for each cell
      obj.matfrac = zeros(obj.grid.ndofs,length(obj.materials.db)-1);
      if isfield(data,"materialFractions")
        frac = load(data.materialFractions);
        obj.matfrac = frac;
      end

      if isfield(data,"materialTag")
        if ischar(data.materialTag)
          loc = str2num(data.materialTag);
        else
          loc = data.materialTag;
        end
        frac = 1/length(loc);
        obj.matfrac(:,loc)=frac;
      end
    end

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
      obj.getState().data.cellDefm = zeros(ndofs,1);
      obj.getState().data.stress = zeros(ndofs,1);
      obj.getState().data.sedmrate = zeros(ndofs,1);
      if obj.nonElasticFlag
        obj.getState().data.voidrate = zeros(ndofs,1);
        obj.getState().data.prestress = obj.getCellsProp('preConStress');
      end
    end

    function prepareSystem(obj)
      % Compute the Neighborhood.
      dofs = obj.grid.getDofsFromIJK;
      ActiveCells = dofs ~= 0;

      actDofs = dofs(ActiveCells);
      actIJK = obj.grid.getIJKfromCellID(actDofs);
      cellNeigh = obj.grid.getNeigh(actIJK);

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

      % Computing the initial porosity
      void0 = obj.getCellsProp('voidrate');
      poro = void0./(1+void0);
      obj.poro0 = poro;
      if obj.nonElasticFlag
        obj.getState().data.voidrate = void0;
      end

      % Initial stress state
      initStres = obj.getCellsProp('initialStress',dofs);
      coords = obj.grid.getCoordCenter(actDofs);

      gamma = obj.materials.getFluid().getSpecificWeight();
      spwg = zeros(length(dofs),1);
      for mat=1:obj.nmat
        spwg = spwg + (obj.materials.getMaterial(mat).SedMaterial.getSpecificWeight() ...
          -gamma)*obj.matfrac(dofs,mat);
      end
      obj.getState().data.stress = -(1-poro).*spwg.*(obj.grid.getColumnMaxHeight()-coords(:,3))-initStres;
    end


    function updateSedRate(obj,dt)
      % Finding the sedimentation rate
      t0 = obj.getStateOld().t;
      sedRate = obj.sedimentHistory.getSedimentationMap(t0,dt);
      obj.getState().data.sedmrate = sedRate;
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
      dt = obj.getState().t-obj.getStateOld().t;
      addSed = dt*obj.getState().data.sedmrate;
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
      %   - Update half-transmissibility
      %   - Extend state vectors
      %   - Update VTK mesh

      lnk=obj.grid;
      newcells = sum(map);

      colheightBG = lnk.columnsHeight;
      obj.matfrac(end+1:end+newcells,:) = sed(map,:)/obj.heightControl;

      colNotTop = colheightBG~=lnk.ncells(3);
      colNotTop = colNotTop(map);

      % Update the grid and locate the new dofs
      newlayer = lnk.grow(map,obj.heightControl);
      [idI,idJ,idK] = lnk.getIJKTop;
      actIJK = [idI(map),idJ(map),idK(map)];
      dofs = lnk.getDofsFromIJK(actIJK);

      % Update initial porosity rate
      voidI = obj.getCellsProp('voidrate',dofs);
      poro = voidI./(1+voidI);
      obj.poro0(end+1:end+newcells) = poro;

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

      % Update the states
      obj.domain.state.data.pressure(end+1:end+newcells) = 0.;
      obj.domain.state.data.cellDefm(end+1:end+newcells) = 0.;

      coords = obj.grid.getCoordCenter(dofs);
      gamma = obj.materials.getFluid().getSpecificWeight();      
      spwg = zeros(length(dofs),1);
      for mat=1:obj.nmat
        spwg = spwg + (obj.materials.getMaterial(mat).SedMaterial.getSpecificWeight() ...
          -gamma)*obj.matfrac(dofs,mat);
      end
      initStres = obj.getCellsProp('initialStress',dofs);
      obj.getState().data.stress(dofs) = -(1-poro).*spwg.*(obj.grid.getColumnMaxHeight(dofs)-coords(dofs,3))-initStres;

      if obj.nonElasticFlag
        obj.domain.state.data.voidrate(end+1:end+newcells) = voidI;
        obj.domain.state.data.prestress(end+1:end+newcells) = obj.getCellsProp('preConStress',dofs);
      end

      % Update the mesh output
      obj.updateMeshOutputOK(map,newlayer);
    end


    function updateMeshOutput(obj,map,newlayer)
      newcells = sum(map);
      if newlayer
        [XX, YY, ZZ] = ndgrid(obj.grid.coordX, obj.grid.coordY, obj.grid.coordZ(end));
        XX = repelem(XX,2,2);
        XX = XX(2:end-1,2:end-1);
        YY = repelem(YY,2,2);
        YY = YY(2:end-1,2:end-1);
        ZZ = repelem(ZZ,2,2);
        ZZ = ZZ(2:end-1,2:end-1);
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

    function updateMeshOutputOK(obj,map,newlayer)
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


    function updateStateElastic(obj,sNew,sOld)
      % Update the Mesh Deformation
      L = obj.grid.getCellHeight();
      L = L(:)+obj.getStateOld().data.cellDefm;
      elasticPerCell = obj.getCellsProp('youngModulus');
      poissonPerCell = obj.getCellsProp('poissonRatio');
      alpha = (1-2*poissonPerCell).*(1+poissonPerCell)./((1-poissonPerCell).*elasticPerCell);
      obj.getState().data.cellDefm=obj.getState().data.cellDefm + alpha.*L.*(sNew-sOld);
    end

    function updateStateNonElastic(obj,sNew,sOld)
      % For the terzaghi column, the deformation is only in the reduction
      % in the void spaces. So the strain is given by
      %          strain = dVoid/(1-void0)

      % Update the Void Rate
      sPre= obj.getState().data.prestress;
      Cc = obj.getCellsProp('compressIdx');
      Cr = obj.getCellsProp('recompressIdx');
      de = SedMaterial.getVoidRatio(abs(sNew),abs(sOld),abs(sPre),Cc,Cr);
      e0 = obj.getStateOld().data.voidrate;
      obj.getState().data.voidrate = e0 + de;

      % Update the Mesh Deformation - vertical deformation (compaction has negative sign)
      L = obj.grid.getCellHeight();
      L = L(:)+obj.getStateOld().data.cellDefm;
      eps = de./(1+e0);      
      obj.getState().data.cellDefm = obj.getStateOld().data.cellDefm + eps.*L;
    end

    function computeHalfTrans(obj)
      % COMPUTEHALFTRANS Computes cells half-transmissibility.
      dofs = (1:obj.grid.ndofs)';
      obj.halfTrans = zeros(obj.grid.ndofs,3);

      % Get mesh dimension and update with the deformation      
      [dx,dy,dz] = obj.grid.getCellsDims();
      dz = dz + obj.getState('cellDefm');

      % Computing half transmissibilities
      condCell = obj.getCellsProp('conductivity');      
      obj.halfTrans(:,1)= (dy.*dz)./(dx/2).*condCell(dofs,1);
      obj.halfTrans(:,2)= (dx.*dz)./(dy/2).*condCell(dofs,2);
      obj.halfTrans(:,3)= (dx.*dy)./(dz/2).*condCell(dofs,3);
    end

    function oedoComp = computeOedometricCompressibility(obj,dofs)
      if ~exist("dofs","var")
        dofs = obj.grid.getActiveDofs;
      end

      sNew = obj.getState().data.stress(dofs);
      if obj.nonElasticFlag
        % Non elastic formulation
        Cc = obj.getCellsProp('compressIdx',dofs);
        void = obj.getState().data.voidrate(dofs);
        oedoComp = -(1/log(10))*(Cc./((1+void).*sNew));
      else
        % Elastic formulation
        elasticPerCell = obj.getCellsProp('youngModulus');
        elasticPerCell = elasticPerCell(dofs);
        poissonPerCell = obj.getCellsProp('poissonRatio');
        poissonPerCell = poissonPerCell(dofs);
        alpha = (1-2*poissonPerCell).*(1+poissonPerCell)./((1-poissonPerCell).*elasticPerCell);
        oedoComp = alpha./(1-alpha.*sNew);
      end      
    end

    function out = getCellsProp(obj,type,dofs)
      if ~exist("dofs","var")
        dofs = 1:obj.grid.ndofs;
      end

      type = lower(type);
      switch type
        case 'conductivity'
          out = zeros(length(dofs),3);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).SedMaterial.getConductivity();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'specificweight'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).SedMaterial.getSpecificWeight();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'voidrate'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).SedMaterial.getVoidRate();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'compressidx'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).SedMaterial.getCompressibilityIdx();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'recompressidx'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).SedMaterial.getReCompressibilityIdx();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'preconstress'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).SedMaterial.getPreConsolidadeStress();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'initialstress'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).SedMaterial.getInitialStress();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'youngmodulus'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).ConstLaw.E;
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'poissonratio'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.materials.getMaterial(mat).ConstLaw.nu;
            out = out + obj.matfrac(dofs,mat).*tmpMat;
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

  end

end