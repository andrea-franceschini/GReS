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
    bcs                                 % Boundary condition
    sedimentHistory = SedimentsMap()    % Time-dependent sediment control

    heightControl double = 0.1          % Height threshold for new cell

    facesNeigh (:,2)                    % Neighboring cell pairs per face
    facesNeighDir (:,1)                 % Face normal direction (1=x,2=y,3=z)
    halfTrans (:,3)                     % Half transmissibility per cell

    H (:,:)
    P (:,:)

    void0 (:,1) double                  % Initial porosity for each cell
    matfrac (:,:)                       % Material fractions per cell
  end

  properties (Access = protected)
    fieldId
  end

  properties (Access = private)
    nmat
    tol = 1e-8;
    niterMx = 100;
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

      % Initialize sediment control
      obj.sedimentHistory = SedimentsMap(input.SedimentMap,obj.nmat,...
        obj.grid.ncells(1:2));

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
    end

    function initialize(obj)
      % Initialize the Sedimentation solver before the simulation starts

      % Set mesh output.
      [coord,conec] = obj.grid.getMesh;
      npts = size(coord,1);
      nelm = size(conec,1);
      obj.mesh.nDim = 3;
      obj.mesh.nCells = nelm;
      obj.mesh.nNodes = npts;
      obj.mesh.cellVTKType = 12*ones(nelm,1);
      obj.mesh.cellNumVerts = 8*ones(nelm,1);
      obj.mesh.cellTag = 8*ones(nelm,1);
      obj.mesh.coordinates = coord;
      obj.mesh.cells = conec;

      % Prepare the faces connectivity
      dofs = (1:nelm)';
      cellNeigh = obj.grid.getNeigh();
      cellsConc = zeros([6*nelm,2]);
      cellsConcDir = ones([6*nelm,1]);
      cellsConcActive = false([6*nelm,1]);
      for ref=1:6
        % Add the connection between the reference cell and the neighborhoods
        tmp = cellNeigh(:,ref);

        % Evaluate the neighbor is smaller than the dof.
        flag = tmp > dofs;
        nfaces = sum(flag);

        cellsConc(ref*nelm+1:ref*nelm+nfaces,:)=[dofs(flag),tmp(flag)];
        cellsConcDir(ref*nelm+1:ref*nelm+nfaces)=ceil(ref/2);
        cellsConcActive(ref*nelm+1:ref*nelm+nfaces)=true;
      end
      obj.facesNeigh = cellsConc(cellsConcActive,:);
      obj.facesNeighDir = cellsConcDir(cellsConcActive);

      % Initialize the States
      obj.getState().data.newcells =0;
      obj.getState().data.pressure = zeros(nelm,1);
      obj.getState().data.cellDefm = zeros(nelm,1);
      obj.getState().data.sedmrate = zeros(nelm,1);
      obj.getState().data.strain = zeros(nelm,1);
      obj.getState().data.stressCons = -abs(obj.getCellsProp('preConStress'));
      obj.getState().data.maxDofUnchanged = 0;
      obj.getState().data.voidrate = zeros(nelm,1);
      obj.getState().data.stress = zeros(nelm,1);
      obj.getState().data.iterTimeStep = 1;
      obj.getState().data.sedimentAcc = zeros(prod(obj.grid.ncells(1:2)),obj.nmat);
      obj.domain.stateOld = copy(obj.domain.getState());

      % Initialize for each cell: porosity, initial stress
      obj.void0 = zeros(nelm,1);
      stStore = zeros(obj.grid.ncells(1:2));
      for nlay = 0:obj.grid.ncells(3)-1
        ddf = obj.grid.dof(:,:,end-nlay);
        map = ddf~=0;
        dtmp = ddf(map);
        fracMat = obj.matfrac(dtmp,:);
        [~,~,dh] = obj.grid.getCellsDims(dtmp);
        stsCons = stStore(map);
        [void,stress] = obj.initializeCell(fracMat,dh/2,stsCons);
        stStore(map) = 2*stress-stsCons;

        obj.void0(dtmp) = void;
        sInit = obj.getCellsProp('initialStress',dtmp);

        obj.getState().data.stress(dtmp) = stress-sInit;
        obj.getStateOld().data.stress(dtmp) = -sInit;
        obj.getState().data.voidrate(dtmp) = void;
        obj.getStateOld().data.voidrate(dtmp) = void;
      end
    end

    function timeStepSetup(obj)
      % initialize the physics solver for the time step

      % Finding the sedimentation rate
      t0 = obj.getStateOld().t;
      dt = obj.domain.state.t-t0;
      sedRate = obj.sedimentHistory.getSedimentationMap(t0,dt);
      obj.getState().data.sedmrate = sedRate;
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
        dofs = obj.grid.getBord(bcId.surface);
        p = getState(obj,"pressure");
        fcorr = false;
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
            fcorr = true;
          otherwise
        end

        switch lower(bcId.type)
          case {'dirichlet'}
            % Correction to the zero pressure to be at the top of the
            % accumulated sediment.
            Trans = obj.halfTrans(dofs,axis);
            if fcorr
              [dx,dy,~] = obj.grid.getCellsDims(dofs);
              sedm=obj.domain.state.data.sedimentAcc;
              dz = sum(sedm,2);
              map = dz~=0;
              if any(map)
                mfrac = sedm(map,:)./dz(map);
                condCell = zeros(sum(map),3);
                for mat=1:obj.nmat
                  tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getConductivity();
                  condCell = condCell + mfrac(:,mat).*tmpMat;
                end
                Tvirt = (dx(map).*dy(map))./(dz(map)).*condCell(:,3);
                Trans(map) = 1./(1./Tvirt+1./Trans(map));
              end
            end

            mu = obj.domain.materials.getFluid().getSpecificWeight();
            dirJ = 1/mu*Trans;
            potential = p(dofs) - bcId.value;
            q = dirJ.*potential;

            nDoF = obj.grid.ndofs;
            bcDofsJ = nDoF*(dofs-1) + dofs;

            obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) = ...
              obj.domain.J{obj.fieldId,obj.fieldId}(bcDofsJ) + dirJ;

            obj.domain.rhs{obj.fieldId}(dofs) = ...
              obj.domain.rhs{obj.fieldId}(dofs) + q;

          case {'neumann'}
            v = bcId.value.*ones(length(dofs),1);
            obj.domain.rhs{obj.fieldId}(dofs) = ...
              obj.domain.rhs{obj.fieldId}(dofs) - sum(vecN.*v,2);

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
      [cellGrow, cellSed] = obj.updateSedAccumulated;

      % Update state for next step
      obj.domain.stateOld = copy(obj.domain.state);
      map = obj.domain.state.data.stress < obj.domain.state.data.stressCons;
      obj.domain.state.data.stressCons(map) = obj.domain.state.data.stress(map);

      % Grow mesh if height threshold is reached
      flagGrow = any(cellGrow);
      if flagGrow
        meshUpdate(obj,cellGrow,cellSed);
      end

      % Update the dof used for the parallel solver
      if flagGrow
        obj.getState().data.maxDofUnchanged = obj.grid.getMaxDofUnchanged;
      end

      % Update the necessity to compute the sedimentation rate.
      obj.getState().data.iterTimeStep = 1;
      obj.getState().data.newcells = sum(cellGrow);
    end

    function updateState(obj,solution)
      % UPDATESTATE Iterate in the time step.

      % Update pressure
      state = getState(obj);
      stateOld = getStateOld(obj);
      state.data.pressure = state.data.pressure + solution;
      dp = state.data.pressure - stateOld.data.pressure;

      % Update the stress state
      dt = state.t-stateOld.t;
      sPrev = stateOld.data.stress;

      dStsMap = obj.getCellsProp('TotalStressVariation');
      dofList = obj.grid.getMapFromDofs();
      dStsMap = dStsMap(dofList);
      if state.data.iterTimeStep == 1
        sCurr = state.data.stress + dp - dt*dStsMap;
      else
        sCurr = sPrev + dp - dt*dStsMap;
      end
      state.data.stress = sCurr;

      % Update the Void Rate
      sCon= state.data.stressCons;
      Cc = obj.getCellsProp('compressIdx');
      Cr = obj.getCellsProp('recompressIdx');
      delta_e = SedimentMaterial.getDeltaVoidRatio(sCurr,sPrev,sCon,Cc,Cr);

      e_prev = stateOld.data.voidrate;
      state.data.voidrate = e_prev + delta_e;

      % Update the Mesh Deformation - vertical deformation
      [~,~,dz] = obj.grid.getCellsDims();
      % dz = dz + obj.getState('cellDefm');  % to compute as eulerian grid
      eps = delta_e./(1+obj.void0);
      strain = stateOld.data.strain + eps;  % <-- Lagrangian strain
      % strain = obj.getStateOld().data.strain + eps./(1+eps); % <-- Eulerian strain

      state.data.strain = strain;
      state.data.cellDefm = strain.*dz;

      % Update the iterator for the time step
      state.data.iterTimeStep = state.data.iterTimeStep + 1;
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
      rhs = rhsStiff + rhsCap;

      K = obj.getCellsProp('TotalStressVariation');
      K = K(obj.grid.getMapFromDofs());
      [dx,dy,dz] = obj.grid.getCellsDims();

      % sedm=obj.domain.state.data.sedimentAcc;
      % vv = dz + obj.getState('cellDefm');
      % tt=obj.grid.getTopDofs();
      % vv(tt) = vv(tt)+sum(sedm,2);
      % volCell = dx.*dy.*vv;
      
      volCell = dx.*dy.*(dz + obj.getState('cellDefm'));
      cb = obj.computeOedometricCompressibility();

      rhsSedm = volCell .* cb .* K;
      rhs = rhs - rhsSedm;
    end

    function computeStiffMat(obj)
      % Compute the Stiffness Matrix
      ncells = obj.grid.ndofs;
      lw = 1/obj.domain.materials.getFluid().getSpecificWeight();

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
      dofs = 1:obj.grid.ndofs;

      % Computing the volumes
      [dx,dy,dz] = obj.grid.getCellsDims();
      volCell = dx.*dy.*(dz + obj.getState('cellDefm'));
      % volCell = dx.*dy.*dz;

      % Computing the storage coefficient
      poro = obj.getState().data.voidrate;
      poro = poro./(1+poro);

      beta = obj.domain.materials.getFluid().getFluidCompressibility();
      oedoComp = obj.computeOedometricCompressibility(dofs);
      PVal = (oedoComp+beta*poro).*volCell;

      % CreatcellDefming the P matrix.
      obj.P = PVal.*speye(obj.grid.ndofs);
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      outPrint = obj.finalizeState(fac,t);
      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function writeSolution(obj,fac,tID)
      outPrint = obj.finalizeState(fac);
      obj.domain.outstate.results(tID).time = outPrint.time;
      obj.domain.outstate.results(tID).pressure = outPrint.pres;
      obj.domain.outstate.results(tID).porosity = outPrint.poro;
      obj.domain.outstate.results(tID).stress = outPrint.stress;
      obj.domain.outstate.results(tID).strain = outPrint.strain;
      obj.domain.outstate.results(tID).compaction = outPrint.comp;
      obj.domain.outstate.results(tID).height = outPrint.height;
    end

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      pointStr = [];

      cellStr = repmat(struct('name', 1, 'data', 1), 2, 1);
      cellStr(1).name = 'pressure';
      cellStr(1).data = state.pres;
      cellStr(2).name = 'piezometric head';
      cellStr(2).data = state.head;

      cellStr(3).name = 'conductivity';
      cellStr(3).data = state.cond;
      cellStr(4).name = 'porosity';
      cellStr(4).data = state.poro;
      cellStr(5).name = 'stress(vertical)';
      cellStr(5).data = state.stress;
      cellStr(6).name = 'strain';
      cellStr(6).data = state.strain;

      cellStr(7).name = 'voidRate';
      cellStr(7).data = state.void;

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

    function prepareMaterialFractions(obj,data)
      % Build the material fractions for each cell
      obj.matfrac = zeros(obj.grid.ndofs,length(obj.domain.materials.solid));
      if isfield(data,"materialFractions")
        frac = load(data.materialFractions);
        obj.matfrac = frac;
      end

      if isfield(data,"materialTag")
        loc = data.materialTag;
        frac = 1/length(loc);
        obj.matfrac(:,loc)=frac;
      end
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

      % 1. Initializate cellSed
      cellSed = zeros(prod(obj.grid.ncells(1:2)),obj.nmat);

      % 2. Update accumulation buffer
      dt = obj.getState().t-obj.getStateOld().t;
      sedAdd = dt*obj.getState().data.sedmrate;
      obj.getState().data.sedimentAcc = obj.getState().data.sedimentAcc + sedAdd;

      % 3. Check for growth trigger (Height >= Threshold)
      colSed = sum(obj.getState().data.sedimentAcc,2);
      cellGrow = colSed/obj.heightControl-1 > -obj.tol; % <-- Very important check

      % 4. Handle overflow for grown cells
      if any(cellGrow)
        gresLog().log(1,"Created %i new cells \n",sum(cellGrow))
        dl = colSed-obj.heightControl;
        sed = (dl./sum(sedAdd,2)).*sedAdd;
        cellSed(cellGrow,:) = obj.getState().data.sedimentAcc(cellGrow,:)-sed(cellGrow,:);
        obj.getState().data.sedimentAcc(cellGrow,:) = sed(cellGrow,:);
      end
    end

    function meshUpdate(obj,map,sed)
      % MESHUPDATE Updates grid topology due to sediment growth.
      %
      % Actions:
      %   - Add new cells
      %   - Update face connectivity
      %   - Update VTK mesh
      %   - Update states

      lnk=obj.grid;
      newcells = sum(map);
      dofs = (obj.grid.ndofs+1:obj.grid.ndofs+newcells)';

      colheightBG = lnk.columnsHeight;
      fracMat = sed(map,:)/obj.heightControl;
      obj.matfrac(end+1:end+newcells,:) = fracMat;

      colNotTop = colheightBG~=lnk.ncells(3);
      colNotTop = colNotTop(map);

      % Update the grid and locate the new dofs
      newlayer = lnk.grow(map,obj.heightControl);

      % Creating the connectives between new cells
      cellNeigh = obj.grid.getNeigh(dofs);

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

      % Update the mesh output
      obj.updateMeshOutput(dofs);

      % Initialize for each cell: porosity, initial stress
      dh = obj.heightControl/2*ones(newcells,1);
      stsCons = zeros(newcells,1);
      [void,stress] = obj.initializeCell(fracMat,dh,stsCons);
      obj.void0(end+1:end+newcells) = void;

      % Update the states
      sInit = obj.getCellsProp('initialStress',dofs);
      obj.domain.state.data.pressure(end+1:end+newcells) = 0.;
      obj.domain.state.data.cellDefm(end+1:end+newcells) = 0.;
      obj.domain.state.data.stressCons(end+1:end+newcells) = -abs(obj.getCellsProp('preConStress',dofs));
      obj.domain.state.data.strain(end+1:end+newcells) = 0.;

      obj.domain.state.data.stress(end+1:end+newcells) = stress-sInit;
      obj.domain.state.data.voidrate(end+1:end+newcells) = void;

      obj.domain.stateOld = copy(obj.domain.getState());
      obj.getStateOld().data.stress(dofs) = -sInit;
    end


    function updateMeshOutput(obj,dofs)
      [coord, conec] = obj.grid.getMesh(dofs);
      npts = size(coord,1);
      nelm = size(conec,1);
      obj.mesh.nCells = obj.mesh.nCells+nelm;
      obj.mesh.nNodes = obj.mesh.nNodes+npts;

      obj.mesh.cellVTKType(end+1:end+nelm) = 12;
      obj.mesh.cellNumVerts(end+1:end+nelm) = 8;
      obj.mesh.cellTag(end+1:end+nelm) = 8;

      obj.mesh.coordinates(end+1:obj.mesh.nNodes,:) = coord;
      obj.mesh.cells(end+1:end+nelm,:) = conec;
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
        dofs = (1:obj.grid.ndofs)';
      end

      % sNew = obj.getState().data.stress(dofs);
      Cc = obj.getCellsProp('compressIdx',dofs);
      Cr = obj.getCellsProp('recompressIdx',dofs);
      sCurr = obj.getState().data.stress(dofs);
      sPrev = obj.getStateOld().data.stress(dofs);
      sCon  = obj.getState().data.stressCons(dofs);
      void = obj.getState().data.voidrate(dofs);
      oedoComp = (1./(1+void)).*SedimentMaterial.getDevVoidRatio(sCurr,sPrev,sCon,Cc,Cr);
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
          dofs = obj.grid.getTopDofs;
          out = zeros(length(dofs),1);
          gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
          omega = obj.getState().data.sedmrate;
          for mat=1:obj.nmat
            gamma_s = obj.domain.materials.getMaterial(mat).ConstLaw.getSpecificWeight();
            out = out + (gamma_s - gamma_w)*omega(:,mat);
          end
          out = (1./(1 + obj.void0(dofs))) .*out;
        otherwise
          out = [];
      end
    end

    function data = finalizeState(obj,fac,t)
      % append state variable to output structure
      sOld = obj.getStateOld();
      sNew = obj.getState();
      if ~exist("t","var")
        t = sNew.t*fac+sOld.t*(1-fac);
      end

      comp = sNew.data.cellDefm*fac+sOld.data.cellDefm*(1-fac);
      data.dl = comp;

      comp = obj.grid.cell2NodeAccByColumnFromBot2Top(comp);
      gamma = obj.domain.materials.getFluid().getSpecificWeight();
      [~,~,zCoordCM] = obj.grid.getCoordCenter();
      voidR = sNew.data.voidrate*fac+sOld.data.voidrate*(1-fac);
      voidR = voidR./(1+voidR);
      % poroA = sNew.data.voidrate./(1+sNew.data.voidrate);
      % poroB = sOld.data.voidrate./(1+sOld.data.voidrate);
      % voidR = fac*poroA+(1-fac)*poroB;

      % mob = (1/obj.domain.materials.getFluid().getDynViscosity());
      % states.flux = computeFlux(obj,p,mob,t);

      data.pres = sNew.data.pressure*fac+sOld.data.pressure*(1-fac);
      data.stress = sNew.data.stress*fac+sOld.data.stress*(1-fac);
      data.strain = sNew.data.strain*fac+sOld.data.strain*(1-fac);
      data.void = sNew.data.voidrate*fac+sOld.data.voidrate*(1-fac);

      data.cond = obj.getCellsProp('conductivity');
      data.comp = comp;
      data.head = zCoordCM+data.pres/gamma;
      data.poro = voidR;
      data.height = obj.grid.getCoordBottom();
      data.time = t;
      % data.height = obj.grid.getCoordTop();
    end

    function [void,stress] = initializeCell(obj,fracMat,dh,stsCons)
      % Computing the e0, gamma_s, Cr, stsInit
      ncells = size(fracMat,1);
      void = zeros(ncells,1);
      gamma_s = zeros(ncells,1);
      gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
      Cr = zeros(ncells,1);
      stsInit = zeros(ncells,1);
      for mat=1:obj.nmat
        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getVoidRate();
        void = void + fracMat(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getSpecificWeight();
        gamma_s = gamma_s + fracMat(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getReCompressibilityIdx();
        Cr = Cr + fracMat(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getInitialStress();
        stsInit = stsInit + fracMat(:,mat).*tmpMat;
      end

      % iterate to find the true value of e and stress
      stsParcial = (gamma_s-gamma_w).*dh;
      % stress = (1./(1+void)).*stsParcial+stsCons;
      [void,stress] = Sedimentation.initialCellProp(void,-stsParcial,...
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