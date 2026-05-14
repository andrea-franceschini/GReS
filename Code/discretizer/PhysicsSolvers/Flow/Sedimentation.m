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

    % void0 (:,1)                         % Initial porosity for each cell
    matfrac (:,:)                       % Material fractions per cell

    deltaStress (:,1)
    voidTop
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
      state = getState(obj);
      ncells = [length(obj.coordX)-1 length(obj.coordY)-1 length(obj.coordZ)-1];

      actDofs = getActiveDof(obj.mdof)';
      state.stressCons(actDofs) = -abs(getCellsProp(obj,'preConStress'));
      stateOld = state;

      % Initialize for each cell
      gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
      stStore = zeros(ncells(1:2));
      for lay = ncells(3):-1:1
        [dof,map] = getDofFromLay(obj.mdof,lay);
        if ~isempty(dof)
          matFract = obj.matfrac(dof,:);
          dh = obj.cellDims(dof,3);
          stsCons = stStore(map);

          ncols = length(dof);
          DStress = zeros(ncols,1);
          voidTmp = zeros(ncols,1);
          stsCellAcc = zeros(ncols,1);
          Cr = zeros(ncols,1);          
          for mat=1:obj.nmat
            gamma_s = obj.domain.materials.getMaterial(mat).ConstLaw.getSpecificWeight();
            DStress = DStress + (gamma_s - gamma_w)*matFract(:,mat).*dh;

            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getVoidRate();
            voidTmp = voidTmp + matFract(:,mat).*tmpMat;

            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getReCompressibilityIdx();
            Cr = Cr + matFract(:,mat).*tmpMat;

            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getInitialStress();
            stsCellAcc = stsCellAcc + matFract(:,mat).*tmpMat;
          end

          % Iterative procedure
          strial = (1./(1+voidTmp)).*DStress+stsCons(map);
          [void,~] = Sedimentation.initialCellProp(voidTmp,-strial,...
            -stsCellAcc,-stsCons(map),Cr,obj.tol,obj.niterMx);
          stress = (1./(1+void)).*DStress+stsCons(map);
          stStore(map) = stress;

          state.stress(dof) = -stress;
          state.voidrate(dof) = void;
          stateOld.stress(dof) = -stress;
          stateOld.voidrate(dof) = void;
        end
      end

      % Set States
      obj.setState(state);
      obj.setStateInit(stateOld);
      obj.setStateOld(stateOld);
    end

    function timeStepSetup(obj)
      % initialize the physics solver for the time step
      dofs = getActiveDof(obj.mdof);
      state = getState(obj);
      stateOld = getStateOld(obj);

      % Finding the sedimentation rate
      t0 = stateOld.time;
      dt = state.time-t0;
      state.sedmRate = getSedimentationMap(obj.sedimentHistory,t0,dt);
      % setState(obj,sedmRate,'sedmRate');

      % Sediments stress contribution
      ncols = prod(obj.mdof.ncells(1:2));      
      matFract = normalize(state.sedmRate,2,'range');
      map = sum(matFract,2)~=0;

      obj.deltaStress = zeros(ncols,1);
      obj.voidTop = zeros(ncols,1);
      Cr = zeros(ncols,1);
      stsInit = zeros(ncols,1);
      gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
      for mat=1:obj.nmat
        gamma_s = obj.domain.materials.getMaterial(mat).ConstLaw.getSpecificWeight();
        obj.deltaStress = obj.deltaStress + (gamma_s - gamma_w)*state.sedmRate(:,mat);

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getVoidRate();
        obj.voidTop = obj.voidTop + matFract(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getReCompressibilityIdx();
        Cr = Cr + matFract(:,mat).*tmpMat;

        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getInitialStress();
        stsInit = stsInit + matFract(:,mat).*tmpMat;
      end

      % iterate to find the true value of e and stress
      stsCons = zeros(ncols,1);
      strial = (1./(1+obj.voidTop(map))).*obj.deltaStress(map)+stsCons(map);
      [void,~] = Sedimentation.initialCellProp(obj.voidTop(map),-strial,...
        -stsInit(map),stsCons(map),Cr(map),obj.tol,obj.niterMx);

      obj.deltaStress(map) = (1./(1+void)).*obj.deltaStress(map,1);
      obj.voidTop(map) = void;

      % state.stress
      dstress = obj.deltaStress(obj.mdof.getMapFromDofs(dofs));
      state.stress(dofs) = state.stress(dofs) - dstress;

      setState(obj,state);
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

      % Get state and set the initial conditions
      state = getState(obj);
      state.iterTimeStep = 1;
      state.strain(:) = 0;
      mapStressCons = state.stress < state.stressCons;
      state.stressCons(mapStressCons) = state.stress(mapStressCons);

      % fprintf('void ratio\n');
      % fprintf('%d\n',state.voidrate(4));
      fprintf('t=%d --- p=%d --- vr=%d\n',state.time,state.pressure(4),state.voidrate(4));
      % fprintf('vr(%d)-%d  --- p(%d)-%d\n',state.time,state.voidrate(4),state.time,state.pressure(4));

      % ------ Update the Sedimentation ------
      % Update the sedimentation 
      dt = state.time-getStateOld(obj,'time');
      sedAdd = dt*state.sedmRate;
      state.sedmAcc = state.sedmAcc + sedAdd;

      newCellSed = zeros(size(state.sedmAcc));
      varCellSed = sedAdd;

      % Check for growth trigger (Height >= Threshold)
      dzSed = sum(state.sedmAcc,2);
      mapNewCells = dzSed > obj.heightControl*(1+obj.minCellHeight); % <-- Very important check

      % Activated the inactive cells
      activeCols = and((dzSed > obj.heightControl*obj.minCellHeight),~state.cellVarAct);
      activetedDofs(obj.mdof,activeCols);
      state.cellVarAct(activeCols)=true;

      % Handle overflow for grown cells
      if any(mapNewCells)
        gresLog().log(1,"Created %i new cells \n",sum(mapNewCells))

        dl = dzSed-obj.heightControl;
        sed = dl.*(sedAdd./sum(sedAdd,2));
        newCellSed(mapNewCells,:) = sed(mapNewCells,:);
        varCellSed(mapNewCells,:) = varCellSed(mapNewCells,:) - newCellSed(mapNewCells,:);
        state.sedmAcc(mapNewCells,:) = sed(mapNewCells,:);
      end

      % ------ Update the cells at the top (before the grow) ------
      dlSedm = sum(sedAdd,2);
      
      dofs = getVarHeightDofs(obj.mdof,state.cellVarAct);
      dzOld = obj.cellDims(dofs,3);
      dlInc = sum(varCellSed(state.cellVarAct,:),2);
      dzNew = dzOld+dlInc;

      fracVarCell = dzOld.*obj.matfrac(dofs,:) + varCellSed(state.cellVarAct,:);
      obj.cellDims(dofs,3) = dzNew;
      obj.matfrac(dofs,:) = fracVarCell./sum(fracVarCell,2);

      state.voidrate(dofs) = (dzOld.*state.voidrate(dofs) + ...
        dlInc.*obj.voidTop(state.cellVarAct))./dzNew;
      state.stress(dofs) = state.stress(dofs) - ...
        obj.deltaStress(state.cellVarAct).*dlInc./dlSedm(state.cellVarAct);

      % ------ Grow the grid if the threshold is reached ------
      newcells = sum(mapNewCells);
      state.newcells = newcells;
      if newcells>0        
        ndofs = obj.mdof.ndofs;
        dofs = (ndofs+1:ndofs+newcells)';
        % obj.matfrac(end+1:end+newcells,:) = normalize(sedNewCell(mapNewCells,:),2,'range');
        obj.matfrac(end+1:end+newcells,:) = newCellSed(mapNewCells,:)./sum(newCellSed(mapNewCells,:),2);

        % Update the grid and locate the new dofs
        newlayer = grow(obj.mdof,mapNewCells);
        if newlayer
          obj.coordZ(end+1) = obj.coordZ(end)+obj.heightControl;
        end

        % Position of the grow
        nlaysByCol = obj.mdof.laysByCol(:);
        atTop = nlaysByCol==obj.mdof.ncells(3);
        newCellNotTop = ~atTop(mapNewCells);

        % Update the cell mesh.
        segmX = diff(obj.coordX);
        segmY = diff(obj.coordY);
        dx = repmat(segmX,[length(segmY) 1]);
        dy = repelem(segmY,length(segmX));
        obj.cellDims(end+1:end+newcells,:) = [dx(mapNewCells) dy(mapNewCells) sum(newCellSed(mapNewCells,:),2)];

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

        dlNew = sum(newCellSed(mapNewCells,:),2);

        % Update the states
        state.pressure(end+1:end+newcells) = 0.;
        state.stress(end+1:end+newcells) = 0.;
        state.strain(end+1:end+newcells) = 0.;
        state.cellDefm(end+1:end+newcells) = 0.;
        state.stressCons(end+1:end+newcells) = 0.;
        state.voidrate(end+1:end+newcells) = 0.;
        state.cellVarAct(mapNewCells) = false;
        state.sedmAcc(mapNewCells,:) = newCellSed(mapNewCells,:);
        state.cellVarAct(mapNewCells) = dlNew > obj.heightControl*obj.minCellHeight;

        % Activate inactive cells
        activetedDofs(obj.mdof,state.cellVarAct);        

        state.voidrate(dofs) = obj.voidTop(mapNewCells);
        state.stress(dofs) = - obj.deltaStress(mapNewCells).*dlNew./dlSedm(mapNewCells);

        state.maxDofUnchanged = getMaxDofUnchanged(obj.mdof);
      end
      
      setStateOld(obj,state);
      setState(obj,state);
    end

    function updateState(obj,solution)
      % UPDATESTATE Iterate in the time step.

      % Update pressure
      dofs = obj.mdof.getActiveDof;
      state = getState(obj);
      stateOld = getStateOld(obj);
      state.pressure(dofs) = state.pressure(dofs) + solution;
      % dp = state.pressure(dofs) - stateOld.pressure(dofs);

      % Update the stress state
      % dt = state.time-stateOld.time;
      sPrev = stateOld.stress(dofs);

      % dStsMap = obj.getCellsProp('TotalStressVariation');
      % dofList = obj.mdof.getMapFromDofs(dofs);
      % dStsMap = dStsMap(dofList);
      % if state.iterTimeStep == 1
      %   sCurr = state.stress(dofs) + dp - dt*dStsMap;
      % else
      %   sCurr = sPrev + dp - dt*dStsMap;
      % end
      % state.stress(dofs) = sCurr;

      sCurr = state.stress(dofs) + solution;
      state.stress(dofs) = sCurr;

      % Update the Void Rate
      sCon = state.stressCons(dofs);
      Cc = obj.getCellsProp('compressIdx',dofs);
      Cr = obj.getCellsProp('recompressIdx',dofs);
      delta_e = SedimentMaterial.getDeltaVoidRatio(sCurr,sPrev,sCon,Cc,Cr);

      e_prev = stateOld.voidrate(dofs);
      state.voidrate(dofs) = e_prev + delta_e;

      % Update the Mesh Deformation - vertical deformation
      % eps = delta_e./(1+obj.void0(dofs));
      eps = delta_e./(1+stateOld.voidrate(dofs));

      % Lagrangian strain model
      % strain = stateOld.strain(dofs) + eps;  
      strain = eps;  
      dz = obj.cellDims(dofs,3)+stateOld.cellDefm(dofs,:);

      % Eulerian strain model
      % strain = stateOld.strain(dofs) + eps./(1+eps);
      % dz = obj.cellDims(dofs,3) + state.cellDefm(dofs);

      state.strain(dofs) = strain;
      state.cellDefm(dofs) = stateOld.cellDefm(dofs) + strain.*dz;

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
      computeStiffMat(obj);
      computeCapMat(obj);
      J = obj.H + obj.P/dt;
    end

    function rhs = computeRhs(obj,dt)
      % Compute the residual of the flow problem

      % get pressure state
      dofs = getActiveDof(obj.mdof);
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());
      p = p(dofs);
      pOld = pOld(dofs);

      rhsStiff = obj.H*p;
      rhsCap = (obj.P/dt)*(p - pOld);
      rhs = rhsStiff + rhsCap;
      
      dl = getStateOld(obj,'cellDefm');
      volCell = prod(obj.cellDims(dofs,1:2),2).*(obj.cellDims(dofs,3) + dl(dofs));
      cb = computeOedometricCompressibility(obj,dofs);

      rhsSedm = volCell .* cb .* obj.deltaStress(getMapFromDofs(obj.mdof,dofs)) ;
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
      dl = getStateOld(obj,'cellDefm');
      volCell = prod(obj.cellDims(dofs,1:2),2).*(obj.cellDims(dofs,3) + dl(dofs));

      % Computing the storage coefficient
      poro = getState(obj,'voidrate');
      poro = poro(dofs)./(1+poro(dofs));

      beta = obj.domain.materials.getFluid().getFluidCompressibility();
      oedoComp = computeOedometricCompressibility(obj,dofs);
      PVal = (oedoComp+beta*poro).*volCell;

      % Create the P matrix.
      obj.P = PVal.*speye(ndofs);
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      % Set mesh output.
      sed = sum(interpolate(obj.domain.state,fac,'sedmAcc'),2);
      comp = interpolate(obj.domain.state,fac,'cellDefm');
      cellVarAct = getState(obj,'cellVarAct');

      obj.grid = makeMeshOutput(obj.mdof,obj.grid,obj.coordX,obj.coordY,...
        obj.coordZ,sed,cellVarAct);
      outPrint = finalizeState(obj,fac);
      outPrint.comp = getComp(obj.mdof,cellVarAct,comp);

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

    function computeHalfTrans(obj)
      % COMPUTEHALFTRANS Computes cells half-transmissibility.
      dofs = getActiveDof(obj.mdof)';
      obj.halfTrans = zeros(obj.mdof.getActiveNdof,3);

      % Get mesh dimension and update with the deformation
      dl = getStateOld(obj,'cellDefm');
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