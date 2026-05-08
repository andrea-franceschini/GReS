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
    mdof                                % Dof Manager for this class
    bcs                                 % Boundary condition
    sedimentHistory = SedimentsMap()    % Time-dependent sediment control

    H (:,:)
    P (:,:)

    facesNeigh (:,2)                    % Neighboring cell pairs per face
    facesNeighDir (:,1)                 % Face normal direction (1=x,2=y,3=z)
    halfTrans (:,3)                     % Half transmissibility per cell
    cellDims (:,3)                      % cell dimensions

    % Grid geometry
    coordX (:,1)                 % X coordinates of grid nodes
    coordY (:,1)                 % Y coordinates of grid nodes
    coordZ (:,1)                 % Z coordinates of grid nodes
    heightControl                % Height threshold for new cell

    void0 (:,1)                         % Initial porosity for each cell
    matfrac (:,:)                       % Material fractions per cell
  end

  properties (Access = protected)
    fieldId
  end

  properties (Access = private)
    nmat
    tol = 1e-8;
    niterMx = 100;
    minCellHeight = 1e-9;
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
      obj.prepareGrid(input.Domain);
      map = [length(obj.coordX) length(obj.coordY) length(obj.coordZ)]-1;

      % Build DOF manager
      obj.mdof = SedDofManager(input.Domain,map);

      % Build the cells dimensions.
      segmX = diff(obj.coordX);
      segmY = diff(obj.coordY);
      segmZ = diff(obj.coordZ);
      [idI,idJ,idK] = obj.mdof.getIJKFromDofs;
      obj.cellDims = [segmX(idI) segmY(idJ) segmZ(idK)];

      % correct the height of the top cell.
      obj.cellDims(obj.mdof.getNotActiveDof,3) = 0.;

      % Build the material fractions for each cell
      obj.prepareMaterialFractions(input.Domain.Initial);

      % Prepare the faces connectivity
      nelm = obj.mdof.ndofs;
      dofs = (1:nelm)';      
      cellNeigh = obj.mdof.getNeigh();
      cellsConc = zeros([6*nelm,2]);
      cellsConcDir = ones([6*nelm,1]);
      cellsConcCreate = false([6*nelm,1]);
      for ref=1:6
        % Add the connection between the reference cell and the neighborhoods
        tmp = cellNeigh(:,ref);

        % Evaluate the neighbor is smaller than the dof.
        flag = tmp > dofs;
        nfaces = sum(flag);

        cellsConc(ref*nelm+1:ref*nelm+nfaces,:)=[dofs(flag),tmp(flag)];
        cellsConcDir(ref*nelm+1:ref*nelm+nfaces)=ceil(ref/2);
        cellsConcCreate(ref*nelm+1:ref*nelm+nfaces)=true;
      end
      obj.facesNeigh = cellsConc(cellsConcCreate,:);
      obj.facesNeighDir = cellsConcDir(cellsConcCreate);

      % Initialize sediment control
      obj.sedimentHistory = SedimentsMap(input.SedimentMap,obj.nmat,map(1:2));

      % Setup BCs
      if ~isfield(input.Boundary,"BC")
        gresLog().error("No boundary was defined for the simulation");
      end
      obj.bcs = input.Boundary.BC(:);

      % Allocate system matrices
      obj.fieldId = 1;
      obj.domain.J{obj.fieldId,obj.fieldId} = [];
      obj.domain.rhs{obj.fieldId} = [];

      % Increase number of variables in the dof manager
      obj.domain.dofm.registerVariable(obj.getField(),entityField.cell,1,1);

      % set the variables for the states.
      state = getState(obj);
      state.newcells = 0;
      state.maxDofUnchanged = 0;
      state.iterTimeStep = 1;

      state.pressure = zeros(nelm,1);
      state.stress = zeros(nelm,1);
      state.strain = zeros(nelm,1);
      state.cellDefm = zeros(nelm,1);
      state.stressCons = zeros(nelm,1);      
      state.voidrate = zeros(nelm,1);

      state.cellVarAct(1:prod(map(1:2)),1) = false;
      state.sedmRate = zeros(prod(map(1:2)),obj.nmat);
      state.sedmAcc  = zeros(prod(map(1:2)),obj.nmat);

      setState(obj,state)
    end

    function initialize(obj)
      % Initialize the Sedimentation solver before the simulation starts

      % Initialize the States
      state = obj.getState;
      ncells = [length(obj.coordX)-1 length(obj.coordY)-1 length(obj.coordZ)-1];

      nelm = obj.mdof.ndofs;      
      actDofs = obj.mdof.getActiveDof';
      state.stressCons(actDofs) = -abs(obj.getCellsProp('preConStress'));
      stateOld = state;

      % Initialize for each cell: porosity, initial stress
      obj.void0 = zeros(nelm,1);
      stStore = zeros(ncells(1:2));
      for lay = ncells(3):-1:1
        [dof,map] = obj.mdof.getDofFromLay(lay);
        if ~isempty(dof)
          fracMat = obj.matfrac(dof,:);
          dh = obj.cellDims(dof,3);
          stsCons = stStore(map);
          [void,stress,~,~] = obj.initializeCell(fracMat,dh/2,stsCons);
          stStore(map) = 2*stress-stsCons;

          obj.void0(dof) = void;
          sInit = obj.getCellsProp('initialStress',dof);

          state.stress(dof) = stress-sInit;
          state.voidrate(dof) = void;

          stateOld.stress(dof) = -sInit;          
          stateOld.voidrate(dof) = void;
        end
      end
      obj.setState(state);
      obj.setStateInit(stateOld);
      obj.setStateOld(stateOld);
    end

    function timeStepSetup(obj)
      % initialize the physics solver for the time step

      % Finding the sedimentation rate
      t0 = obj.getStateOld('time');
      dt = obj.getState('time')-t0;
      sedmRate = obj.sedimentHistory.getSedimentationMap(t0,dt);
      obj.setState(sedmRate,'sedmRate');
    end

    function assembleSystem(obj,dt)
      % Update the cells half transmissibility
      obj.computeHalfTrans();

      % Assembling the system.
      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);
    end

    function applyBC(obj,bcId,t)
      for bcId = obj.bcs
        dofs = obj.mdof.getBord(bcId.surface);
        dofsR = obj.mdof.reScale(dofs)';
        p = getState(obj,"pressure");
        switch lower(bcId.surface)
          case {"xmin","xmax"}
            axis=1;
            vecN = [1 0 0];
          case {"ymin","ymax"}
            axis=2;
            vecN = [0 1 0];
          case {"zmin"}
            axis=3;
            vecN = [0 0 1];
          case {"zmax"}
            axis=3;
            vecN = [0 0 1];
          otherwise
        end

        switch lower(bcId.type)
          case {'dirichlet'}
            trans = obj.halfTrans(dofsR,axis);
            mu = obj.domain.materials.getFluid().getSpecificWeight();
            dirJ = 1/mu*trans;
            potential = p(dofs) - bcId.value;
            q = dirJ.*potential;

            nDoF = obj.mdof.getActiveNdof;
            bcDofsJ = nDoF*(dofsR-1) + dofsR;

            obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) = ...
              obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) + dirJ;

            obj.domain.rhs{obj.fieldId}(dofsR) = ...
              obj.domain.rhs{obj.fieldId}(dofsR) + q;

          case {'neumann'}
            v = bcId.value.*ones(length(dofs),1);
            obj.domain.rhs{obj.fieldId}(dofsR) = ...
              obj.domain.rhs{obj.fieldId}(dofsR) - sum(vecN.*v,2);

          otherwise
            gresLog().error("Error in %s: Boundary condition type '%s' is not " + ...
              "available in applyBC()",class(obj),bcType)
        end
      end
    end

    function applyDirVal(obj,bcId,t)
      return
    end

    function advanceState(obj)
      % ADVANCESTATE Finalizes time step and updates grid topology.

      % Update sediment accumulation
      [cellGrow, cellSed, sedmCols, activeCol] = updateSedAccumulated(obj);

      % Update state for next step
      state = getState(obj);
      setStateOld(obj,state);

      map = state.stress < state.stressCons;
      state.stressCons(map) = state.stress(map);
      setState(obj,state.stressCons,'stressCons');      
      setState(obj,sedmCols,'sedmAcc');

      % Active and Update the cells with variable height
      % if ~any(cellGrow)
      updateTopCells(obj,activeCol,true);
      % end

      % Grow mesh if height threshold is reached
      if any(cellGrow)
        meshUpdate(obj,cellGrow,cellSed);

        % Update the dof used for the parallel solver
        setState(obj,obj.mdof.getMaxDofUnchanged,'maxDofUnchanged');
      end

      % Update the necessity to compute the sedimentation rate.
      setState(obj,1,'iterTimeStep');
      setState(obj,sum(cellGrow),'newcells');
    end

    function updateState(obj,solution)
      % UPDATESTATE Iterate in the time step.

      % Update pressure
      dofs = obj.mdof.getActiveDof;
      state = getState(obj);
      stateOld = getStateOld(obj);
      state.pressure(dofs) = state.pressure(dofs) + solution;
      dp = state.pressure(dofs) - stateOld.pressure(dofs);

      % Update the stress state
      dt = state.time-stateOld.time;
      sPrev = stateOld.stress(dofs);

      dStsMap = obj.getCellsProp('TotalStressVariation');
      dofList = obj.mdof.getMapFromDofs(dofs);
      dStsMap = dStsMap(dofList);
      if state.iterTimeStep == 1
        sCurr = state.stress(dofs) + dp - dt*dStsMap;
      else
        sCurr = sPrev + dp - dt*dStsMap;
      end
      state.stress(dofs) = sCurr;

      % Update the Void Rate
      sCon = state.stressCons(dofs);
      Cc = obj.getCellsProp('compressIdx',dofs);
      Cr = obj.getCellsProp('recompressIdx',dofs);
      delta_e = SedimentMaterial.getDeltaVoidRatio(sCurr,sPrev,sCon,Cc,Cr);

      e_prev = stateOld.voidrate(dofs);
      state.voidrate(dofs) = e_prev + delta_e;

      % Update the Mesh Deformation - vertical deformation
      eps = delta_e./(1+obj.void0(dofs));

      % Lagrangian strain model
      strain = stateOld.strain(dofs) + eps;  
      dz = obj.cellDims(dofs,3);

      % Eulerian strain model
      % strain = stateOld.strain(dofs) + eps./(1+eps);
      % dz = obj.cellDims(dofs,3) + state.cellDefm(dofs);

      state.strain(dofs) = strain;
      state.cellDefm(dofs) = strain.*dz;

      % Update the iterator for the time step
      state.iterTimeStep = state.iterTimeStep + 1;

      obj.setState(state);
    end

    function goBackState(obj,varargin)
      % base method to move back the state when convergence is not reached
      dofs = obj.mdof.getBordZMax;
      sedm = normalize(obj.getStateOld("sedmAcc"), 2,'range');      
      obj.matfrac(dofs,:)=sedm;
      cell = obj.getStateOld("cellVarAct") ~= obj.getState("cellVarAct");

      if any(cell)
        obj.mdof.deactiveDofs(cell);
      end

      obj.domain.state = copy(obj.domain.stateOld);
    end

    function J = computeMat(obj,dt)
      % Recompute elementary matrices
      obj.computeStiffMat;
      obj.computeCapMat;
      J = obj.H + obj.P/dt;
    end

    function rhs = computeRhs(obj,dt)
      % Compute the residual of the flow problem

      % get pressure state
      dofs = obj.mdof.getActiveDof;
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());
      p = p(dofs);
      pOld = pOld(dofs);

      rhsStiff = obj.H*p;
      rhsCap = (obj.P/dt)*(p - pOld);
      rhs = rhsStiff + rhsCap;

      K = obj.getCellsProp('TotalStressVariation');
      K = K(obj.mdof.getMapFromDofs(dofs));
      dl = obj.getState('cellDefm');
      volCell = prod(obj.cellDims(dofs,1:2),2).*(obj.cellDims(dofs,3) + dl(dofs));
      cb = obj.computeOedometricCompressibility(dofs);

      rhsSedm = volCell .* cb .* K;
      rhs = rhs - rhsSedm;
    end

    function computeStiffMat(obj)
      % Compute the Stiffness Matrix
      dofs = obj.mdof.getActiveDof;
      ndofs = obj.mdof.getActiveNdof;
      notdof = obj.mdof.getNotActiveDof;
      lw = 1/obj.domain.materials.getFluid().getSpecificWeight();

      actFaces = ~sum(ismember(obj.facesNeigh,notdof),2);
      [~, neight] = ismember(obj.facesNeigh(actFaces,:),dofs);

      idi = sub2ind([ndofs 3], neight(:,1), obj.facesNeighDir(actFaces));
      idk = sub2ind([ndofs 3], neight(:,2), obj.facesNeighDir(actFaces));
      Tii = obj.halfTrans(idi);
      Tik = obj.halfTrans(idk);
      trans = lw./(1./Tii+1./Tik);

      sumDiagTrans = accumarray(neight(:),repmat(trans,[2,1]),[ndofs,1]);
      obj.H = sparse([neight(:,1); neight(:,2); (1:ndofs)'],...
        [neight(:,2); neight(:,1); (1:ndofs)'],...
        [-trans; -trans; sumDiagTrans], ndofs, ndofs);
    end

    function computeCapMat(obj)
      % COMPUTECAPMAT Builds the storage matrix (P).
      %
      % Includes:
      %   - Rock compressibility
      %   - Fluid compressibility

      % Selecting the active dofs
      dofs = obj.mdof.getActiveDof;
      ndofs = obj.mdof.getActiveNdof;

      % Computing the volumes
      dl = obj.getState('cellDefm');
      volCell = prod(obj.cellDims(dofs,1:2),2).*(obj.cellDims(dofs,3) + dl(dofs));

      % Computing the storage coefficient
      poro = obj.getState('voidrate');
      poro = poro(dofs)./(1+poro(dofs));

      beta = obj.domain.materials.getFluid().getFluidCompressibility();
      oedoComp = obj.computeOedometricCompressibility(dofs);
      PVal = (oedoComp+beta*poro).*volCell;

      % Create the P matrix.
      obj.P = PVal.*speye(ndofs);
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      % Set mesh output.
      sed = sum(interpolate(obj.domain.state,fac,'sedmAcc'),2);
      comp = interpolate(obj.domain.state,fac,'cellDefm');
      cellVarAct = getState(obj,'cellVarAct');

      obj.grid = obj.mdof.makeMeshOutput(obj.grid,obj.coordX,obj.coordY,...
        obj.coordZ,sed,cellVarAct);
      outPrint = obj.finalizeState(fac);
      outPrint.comp = obj.mdof.getComp(cellVarAct,comp);

      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function writeSolution(obj,fac,tID)
      outPrint = finalizeState(obj,fac);
      comp = interpolate(obj.domain.state,fac,'cellDefm');
      cellVarAct = getState(obj,'cellVarAct');

      obj.domain.outstate.results(tID).time = outPrint.time;
      obj.domain.outstate.results(tID).pressure = outPrint.pres;
      obj.domain.outstate.results(tID).porosity = outPrint.poro;
      obj.domain.outstate.results(tID).stress = outPrint.stress;
      obj.domain.outstate.results(tID).strain = outPrint.strain;
      obj.domain.outstate.results(tID).compaction = obj.mdof.getComp(cellVarAct,comp);
    end

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      cellStr = repmat(struct('name', 1, 'data', 1), 6, 1);
      cellStr(1).name = 'pressure';
      cellStr(1).data = state.pres;
      cellStr(2).name = 'stress';
      cellStr(2).data = state.stress;
      cellStr(3).name = 'strain';
      cellStr(3).data = state.strain;
      cellStr(4).name = 'voidRate';
      cellStr(4).data = state.void;
      cellStr(5).name = 'porosity';
      cellStr(5).data = state.poro;
      cellStr(6).name = 'conductivity';
      cellStr(6).data = state.cond;

      pointStr(1).name = 'compaction';
      pointStr(1).data = state.comp;
    end

    function out = isLinear(obj)
      out = false;
    end
  end

  methods (Access = private)

    function checkInput(obj,input)
      flag = false;

      % Check the material is well defined.
      obj.nmat = length(obj.domain.materials.solid);

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

    function prepareGrid(obj,data)
      % Constructing the grid      
      if strcmp(data.Grid.type,"classic")
        ncells = data.Grid.division;
        dim = data.Grid.size;
        obj.coordX = linspace(0,dim(1),ncells(1)+1);
        obj.coordY = linspace(0,dim(2),ncells(2)+1);
        obj.coordZ = linspace(0,dim(3),ncells(3)+1);
      elseif strcmp(data.Grid.type,"explicit")
        % Internal initialization of grid, maps, and material layers.
        if isfield(data,'xfile')
          obj.coordX = load(data.xfile);
        else
          obj.coordX = data.xcoord;
        end

        if isfield(data,'yfile')
          obj.coordY = load(data.yfile);
        else
          obj.coordY = data.ycoord;
        end

        if isfield(data,'zfile')
          obj.coordZ = load(data.zfile);
        else
          obj.coordZ = data.zcoord;
        end
      else
        gresLog().error("The grid parameters for the simulation is not well defined!");
      end
      obj.heightControl = data.NewCellHeightControl;

      obj.coordZ(end+1)=obj.coordZ(end)+obj.heightControl;
    end

    function prepareMaterialFractions(obj,data)
      % Build the material fractions for each cell
      ndofs = obj.mdof.getActiveNdof;
      obj.matfrac = zeros(ndofs,length(obj.domain.materials.solid));
      if isfield(data,"materialFractions")
        frac = load(data.materialFractions);
        obj.matfrac = frac;
      end

      if isfield(data,"materialTag")
        loc = data.materialTag;
        frac = 1/length(loc);
        obj.matfrac(:,loc)=frac;
      end

      % Add the cell the space for the inactive cell.
      inDofs = length(obj.mdof.getNotActiveDof);
      obj.matfrac(end+1:end+inDofs,:)=0.;
    end

    function [cellGrow, cellSed, sedAcc, activeCol] = updateSedAccumulated(obj)
      % UPDATESEDACCUMULATED Updates sediment buffer and triggers grid growth.
      %
      %   Calculates new deposition, updates the accumulation buffer, and
      %   identifies columns reaching the height threshold for new cells.
      %
      % Outputs:
      %   cellGrow - Logical array indicating columns that grow
      %   cellSed  - Sediment assigned to newly created cells
      %   sedmCols - Total amount of sediments in each column
      %  activeCol - The columns of which the cell with variable height is
      %              activeted.

      % 1. Initialization cellSed
      sedAcc = getState(obj,'sedmAcc');
      cellSed = zeros(size(sedAcc));

      % 2. Update accumulation buffer
      dt = getState(obj,'time')-getStateOld(obj,'time');
      sedAdd = dt*getState(obj,'sedmRate');
      sedAcc = sedAcc + sedAdd;

      % 3. Check for growth trigger (Height >= Threshold)
      actVarCell = getState(obj,"cellVarAct");
      colSed = sum(sedAcc,2);
      cellGrow = colSed/obj.heightControl-1 > -obj.tol; % <-- Very important check

      % 4. Check the inactive cells with enough sediment to be activated 
      activeCol = and((colSed > obj.heightControl*obj.minCellHeight),~actVarCell);

      % 5. Handle overflow for grown cells
      if any(cellGrow)
        gresLog().log(1,"Created %i new cells \n",sum(cellGrow))
        dl = colSed-obj.heightControl;
        sed = (dl./sum(sedAdd,2)).*sedAdd;
        cellSed(cellGrow,:) = sed(cellGrow,:);
        sedAcc(cellGrow,:) = sedAcc(cellGrow,:)-sed(cellGrow,:);
      end
    end

    function updateTopCells(obj,activeCol,flag)
      state = getState(obj);
      stateOld = getStateOld(obj);

      sedm = state.sedmAcc;
      sedmOld = stateOld.sedmAcc;
      dsedm = sedm - sedmOld;

      % Activate inactive cells
      activetedDofs(obj.mdof,activeCol);
      state.cellVarAct(activeCol) = true;
      setState(obj,state.cellVarAct,'cellVarAct');

      actCells = state.cellVarAct;
      dofs = getVarHeightDofs(obj.mdof,actCells);

      % Change the height of the tops cells.
      % if flag
      %   % Fraction of each column
      %   fracSOld = sum(sedmOld(actCells,:),2)/sum(sedm(actCells,:),2);
      %   fracDSed = sum(dsedm(actCells,:),2)/sum(sedm(actCells,:),2);
      %   dh = sum(dsedm(actCells,:),2);
      % else
        dh = sum(sedm(actCells,:),2);
      % end
      % dh = sum(sedm(actCells,:),2);
      obj.cellDims(dofs,3) = sum(sedm(actCells,:),2);

      % Initialize for each cell: porosity, initial stress
      stsCons = zeros(size(dh));
      obj.matfrac(dofs,:) = normalize(sedm(actCells,:), 2,'range');
      [void,stress,void0,stress0] = initializeCell(obj,obj.matfrac(dofs,:),dh/2,stsCons);
      % if and(flag,length(obj.void0)>max(dofs))
      %   void = 1./(1./fracSOld.*obj.void0(dofs)+1./fracDSed.*void);
      % end
      obj.void0(dofs) = void;

      % Update the states
      sInit = getCellsProp(obj,'initialStress',dofs);
      state.stressCons(dofs) = -abs(getCellsProp(obj,'preConStress',dofs));
      state.voidrate(dofs) = void;
      % if flag
        % state.stress(dofs) = -stress;
      % else
        state.stress(dofs) = stress-sInit;
      % end
      
      stateOld.voidrate(dofs) = void;
      stateOld.stress(dofs) = -sInit;

      setState(obj,state);
      setStateOld(obj,stateOld);
    end

    function meshUpdate(obj,map,sed)
      % MESHUPDATE Updates grid topology due to sediment growth.
      %
      % Actions:
      %   - Add new cells
      %   - Update face connectivity
      %   - Update VTK mesh
      %   - Update states

      newcells = sum(map);
      ndofs = obj.mdof.ndofs;
      dofs = (ndofs+1:ndofs+newcells)';
      obj.matfrac(end+1:end+newcells,:) = zeros(newcells,obj.nmat);

      % Update the grid and locate the new dofs
      newlayer = grow(obj.mdof,map);
      if newlayer
        obj.coordZ(end+1) = obj.coordZ(end)+obj.heightControl;
      end

      % Position of the grow
      nlaysByCol = obj.mdof.laysByCol(:);
      atTop = nlaysByCol==obj.mdof.ncells(3);
      newCellNotTop = ~atTop(map);

      % Update the cell mesh.
      segmX = diff(obj.coordX);
      segmY = diff(obj.coordY);
      dx = repmat(segmX,[length(segmY) 1]);
      dy = repelem(segmY,length(segmX));
      obj.cellDims(end+1:end+newcells,:) = [dx(map) dy(map) sum(sed(map,:),2)];

      % Creating the connectives between new cells
      neigh = getNeigh(obj.mdof,dofs);

      face = zeros([5*newcells,2]);
      faceAxis = ones([5*newcells,1]);
      faceAct = false([5*newcells,1]);

      face(1:newcells,:)=[neigh(:,5) dofs];
      faceAxis(1:newcells)=3;
      faceAct(1:newcells)=true;
      for ref=1:4
        % Add the connection between the reference cell and the neighborhoods
        tmp = neigh(:,ref);

        % Evaluate the neighbor if the column is not at the top.
        flag0 = and(and(tmp<dofs,tmp~=0),newCellNotTop);

        % Evaluate the neighbor if the column is at the top.
        flag1 = tmp > dofs;

        % Combine the two case.
        flag = or(flag0,flag1);
        nfaces = sum(flag);

        face(ref*newcells+1:ref*newcells+nfaces,:)=[dofs(flag),tmp(flag)];
        faceAxis(ref*newcells+1:ref*newcells+nfaces)=ceil(ref/2);
        faceAct(ref*newcells+1:ref*newcells+nfaces)=true;
      end
      nfaces = sum(faceAct);
      obj.facesNeigh(end+1:end+nfaces,:) = face(faceAct,:);
      obj.facesNeighDir(end+1:end+nfaces,:) = faceAxis(faceAct);

      % Update the states
      state = getState(obj);
      stateOld = getStateOld(obj);
      
      state.pressure(end+1:end+newcells) = 0.;
      state.stress(end+1:end+newcells) = 0.;
      state.strain(end+1:end+newcells) = 0.;
      state.cellDefm(end+1:end+newcells) = 0.;
      state.stressCons(end+1:end+newcells) = 0.;
      state.voidrate(end+1:end+newcells) = 0.;
      state.cellVarAct(map) = false;
      state.sedmAcc(map,:) = sed(map,:);
      setState(obj,state);

      stateOld.pressure(end+1:end+newcells) = 0.;
      stateOld.stress(end+1:end+newcells) = 0.;
      stateOld.strain(end+1:end+newcells) = 0.;
      stateOld.cellDefm(end+1:end+newcells) = 0.;
      stateOld.stressCons(end+1:end+newcells) = 0.;
      stateOld.voidrate(end+1:end+newcells) = 0.;
      stateOld.cellVarAct(map) = false;
      stateOld.sedmAcc(map,:) = sed(map,:);
      setStateOld(obj,stateOld);

      % Activeted the new cell if necessary
      activeCol = state.cellVarAct;
      activeCol(map) = sum(sed(map,:),2) > obj.heightControl*obj.minCellHeight;

      obj.updateTopCells(activeCol,false);
    end

    function computeHalfTrans(obj)
      % COMPUTEHALFTRANS Computes cells half-transmissibility.
      dofs = getActiveDof(obj.mdof)';
      obj.halfTrans = zeros(obj.mdof.getActiveNdof,3);

      % Get mesh dimension and update with the deformation
      dl = getState(obj,'cellDefm');
      dx = obj.cellDims(dofs,1);
      dy = obj.cellDims(dofs,2);
      dz = obj.cellDims(dofs,3) + dl(dofs);

      % Computing half transmissibilities
      condCell = getCellsProp(obj,'conductivity',dofs);
      obj.halfTrans(:,1)= (dy.*dz)./(dx/2).*condCell(:,1);
      obj.halfTrans(:,2)= (dx.*dz)./(dy/2).*condCell(:,2);
      obj.halfTrans(:,3)= (dx.*dy)./(dz/2).*condCell(:,3);
    end

    function oedoComp = computeOedometricCompressibility(obj,dofs)
      state = getState(obj);
      stressOld = getStateOld(obj,'stress');      
      Cc = getCellsProp(obj,'compressIdx',dofs);
      Cr = getCellsProp(obj,'recompressIdx',dofs);
      sCurr = state.stress(dofs);
      sPrev = stressOld(dofs);
      sCon  = state.stressCons(dofs);
      void = state.voidrate(dofs);
      oedoComp = (1./(1+void)).*SedimentMaterial.getDevVoidRatio(sCurr,sPrev,sCon,Cc,Cr);
    end

    function out = getCellsProp(obj,type,dofs)
      if ~exist("dofs","var")
        dofs = getActiveDof(obj.mdof);
      end

      type = lower(type);
      switch type
        case 'conductivity'
          out = zeros(length(dofs),3);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getConductivity();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
          % try to impose kz0*10^((e-e0)/Ck)
          % dvoid = obj.getState("voidrate")-obj.void0;
          % out(:,3)=out(:,3).*10.^(dvoid/10);
          % out(:,1:2)=[3*out(:,3) 3*out(:,3)];
        case 'compressidx'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getCompressibilityIdx();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'recompressidx'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getReCompressibilityIdx();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'preconstress'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getPreConsolidadeStress();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'initialstress'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getInitialStress();
            out = out + obj.matfrac(dofs,mat).*tmpMat;
          end
        case 'totalstressvariation'
          % dofs = obj.grid.getTopDofs;
          dofs = getBordZMax(obj.mdof);
          out = zeros(length(dofs),1);
          gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
          omega = getState(obj,'sedmRate');
          for mat=1:obj.nmat
            gamma_s = obj.domain.materials.getMaterial(mat).ConstLaw.getSpecificWeight();
            out = out + (gamma_s - gamma_w)*omega(:,mat);
          end
          out = (1./(1 + obj.void0(dofs))) .*out;
        otherwise
          out = [];
      end
    end

    function data = finalizeState(obj,fac)
      % append state variable to output structure
      dofs = getActiveDof(obj.mdof)';
      state = interpolate(obj.domain.state,fac);
      voidR = state.voidrate(dofs)./(1+state.voidrate(dofs));

      data.time = state.time;
      data.pres = state.pressure(dofs);
      data.stress = state.stress(dofs);
      data.strain = state.strain(dofs);
      data.void = state.voidrate(dofs);
      data.cond = getCellsProp(obj,'conductivity');
      data.poro = voidR;
    end

    function [void,stress,void0,stress0] = initializeCell(obj,fracMat,dh,stsCons)
      % Computing the e0, gamma_s, Cr, stsInit
      ncells = size(fracMat,1);
      void0 = zeros(ncells,1);
      gamma_s = zeros(ncells,1);
      gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
      Cr = zeros(ncells,1);
      stsInit = zeros(ncells,1);
      for mat=1:obj.nmat
        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getVoidRate();
        void0 = void0 + fracMat(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getSpecificWeight();
        gamma_s = gamma_s + fracMat(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getReCompressibilityIdx();
        Cr = Cr + fracMat(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getInitialStress();
        stsInit = stsInit + fracMat(:,mat).*tmpMat;
      end

      % iterate to find the true value of e and stress
      stsParcial = (gamma_s-gamma_w).*dh;
      stress0 = (1./(1+void0)).*stsParcial+stsCons;
      [void,stress] = Sedimentation.initialCellProp(void0,-stsParcial,...
        -stsInit,stsCons,Cr,obj.tol,obj.niterMx);      
    end

  end

  methods (Static)

    function [void, stress] = initialCellProp(voidInit,stsParcial,stsInit,stsCons,Cr,tolerance,maxNIter)
      % INITIALCELLPROP function to initialize the cell void rate and stress
      % using a iterative procedure.
      %
      % $$ \sigma_{trial} = (\frac{1}{1+e})*\sigma_{parcial} + cte $$,
      %
      % $$ e_{c} = e_{0} - C_{r} \log(\frac{\sigma_{trial}}{\sigma_{init}}) $$
      %
      % where:
      % - $$e_{0}$$ - the initial void rate
      % - $$e_{c}$$ - the current void rate
      % - $$C_{r}$$ - the re compression index
      % - $$\sigma_{init}$$ - the initial stress
      % - $$\sigma_{trial}$$ - the trial stress
      % - $$cte$$ - a constant stress value
      %
      % The iterative process consist in evaluate the stress using the void
      % rate and the void rate using the stress. This procedure is maintain
      % util the difference of void rate is smaller than a tolerance value
      % between iterations.
      error = 1e3*tolerance;
      void = voidInit;
      iter = 1;
      dvoidOld = zeros(length(void),1);
      while and(error > tolerance, iter<maxNIter)
        stress = (1./(1+void)).*stsParcial+stsCons;
        dvoid = -Cr.*log10(stress./stsInit);
        error = max(abs((dvoid-dvoidOld)./void));
        void = voidInit + dvoid;
        dvoidOld = dvoid;
        iter=iter+1;
      end
    end

    function out = getField()
      out = "pressure";
    end

    function out = isSymmetric()
      out = true;
    end

  end

end