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

    cellMatBranch
    cbFreeze = false
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

      % Build Grid/DOF manager
      obj.mdof = domainManagerSed(input.Domain,obj.nmat);
      ncells = obj.mdof.ndofs;
      surfdim = size(obj.mdof.cellsByCol);

      % Initialize sediment control
      obj.sedimentHistory = SedimentsMap(input.SedimentMap,obj.nmat,surfdim);

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

      state.pressure = zeros(ncells,1);
      state.stress = zeros(ncells,1);
      state.strain = zeros(ncells,1);
      state.strainAcc = zeros(ncells,1);
      state.cellDefm = zeros(ncells,1);
      state.Sp = zeros(ncells,1);      
      state.e = zeros(ncells,1);
      state.e0 = zeros(ncells,1);
      state.age = zeros(ncells,1);      
      state.cb = zeros(ncells,1);

      state.sedmRate = zeros(ncells,obj.nmat);
      state.sedmAcc  = zeros(ncells,obj.nmat);

      setState(obj,state)
    end

    function initialize(obj)
      % Initialize the Sedimentation solver before the simulation starts

      % Initialize the States
      state = getState(obj);
      state.age(:) = state.time;      

      % Initialize for each cell
      ncells = obj.mdof.ndofs;
      matTag = find(max(obj.mdof.matfrac,[],1));
      gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
      gamma_s = obj.domain.materials.getMaterial(matTag).ConstLaw.gamma;

      void0 = obj.domain.materials.getMaterial(matTag).ConstLaw.e0*ones(ncells,1);
      Cc = obj.domain.materials.getMaterial(matTag).ConstLaw.Cc*ones(ncells,1);
      Cr = obj.domain.materials.getMaterial(matTag).ConstLaw.Cr*ones(ncells,1);
      Sp = obj.domain.materials.getMaterial(matTag).ConstLaw.Sp*ones(ncells,1);

      % Iterative procedure
      cellWeight = (gamma_s - gamma_w)*diff(obj.mdof.coordZ);
      Strial = cellWeight*(1./(1+void0));
      [void,~] = Sedimentation.initialCellProp(void0,-Strial,-Sp,...
        Cr,obj.tol,obj.niterMx);

      Stress = cellWeight*(1./(1+void));
      mapPreCon = Sp<Stress;
      Sp(mapPreCon)=Stress(mapPreCon);

      state.stress(:) = -Stress;
      state.Sp(:) = -Sp;
      state.e(:) = void;
      state.e0 = SedimentMaterial.getVoidPreCon(Stress,Sp,void,Cr,Cc);

      % Set States
      obj.setState(state);
      obj.setStateInit(state);
      obj.setStateOld(state);
    end

    function timeStepSetup(obj)
      % initialize the physics solver for the time step
      lnk = obj.mdof;
      dofs = 1:lnk.ndofs;      
      state = getState(obj);
      stateOld = getStateOld(obj);      

      % ncols = prod(obj.mdof.ncells(1:2));
      topDim = size(lnk.cellsByCol);
      ncols = numel(lnk.cellsByCol);
      lnk.top_dS(:) = 0.;
      lnk.top_e(:) = 0.;
      lnk.top_e0(:) = 0.;

      % Finding the sedimentation rate
      t0 = stateOld.time;
      dt = state.time-t0;
      state.sedmRate = getSedimentationMap(obj.sedimentHistory,t0,dt);
      sedm = state.sedmRate;
      dlSedm = sum(sedm,2);

      matFract = sedm./dlSedm;
      matFract(dlSedm==0,:) = 0;
      map = sum(matFract,2)~=0;
      
      Cc = zeros(ncols,1);
      Cr = zeros(ncols,1);
      stsInit = zeros(ncols,1);
      sedmWeight = zeros(ncols,1);
      gamma_w = obj.domain.materials.getFluid().getSpecificWeight();
      for mat=1:obj.nmat
        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw;

        sedmWeight = sedmWeight + (tmpMat.gamma - gamma_w)*sedm(:,mat);
        lnk.top_e = lnk.top_e + reshape(matFract(:,mat),topDim).*tmpMat.e0;
        Cc = Cc + matFract(:,mat).*tmpMat.Cc;
        Cr = Cr + matFract(:,mat).*tmpMat.Cr;
        stsInit = stsInit + matFract(:,mat).*tmpMat.S0;
      end

      % iterate to find the true value of e and stress
      strial = (1./(1+lnk.top_e(map))).*sedmWeight(map);
      [void,~] = Sedimentation.initialCellProp(lnk.top_e(map),-strial,...
        -stsInit(map),Cr(map),obj.tol,obj.niterMx);

      lnk.top_dS(map) = (1./(1+void)).*sedmWeight(map);
      lnk.top_e(map) = void;
      lnk.top_e0(map) = SedimentMaterial.getVoidPreCon(lnk.top_dS(map),stsInit(map),void,Cr(map),Cc(map));

      dstress = lnk.top_dS(getMapFromDofs(lnk,dofs))';
      state.stress(dofs) = state.stress(dofs) - dstress*dt;

      setState(obj,state);
    end

    function assembleSystem(obj,dt)
      % Update the cells half transmissibility
      obj.MaterialMap();
      % obj.computeHalfTrans();
      obj.computeCb();

      % Assembling the system.
      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);
    end

    function applyBC(obj,bcId,t)
      for bcId = obj.bcs        
        p = getState(obj,"pressure");
        dofs = getBord(obj.mdof,bcId.surface);
        switch lower(bcId.surface)
          case {"xmin","xmax"}
            vecN = [1 0 0];
          case {"ymin","ymax"}
            vecN = [0 1 0];
          case {"zmin","zmax"}
            vecN = [0 0 1];
          otherwise
            gresLog().error("Surface to apply the boundary condition is not allowed");
        end

        switch lower(bcId.type)
          case {'dirichlet'}
            trans = getBcHalfTrans(obj.mdof,find(vecN),...
              getCellsProp(obj,'conductivity',dofs),...
              getStateOld(obj,'cellDefm'),dofs);
            lw = 1/obj.domain.materials.getFluid().getSpecificWeight();
            dirJ = lw*trans;
            potential = p(dofs) - bcId.value;
            q = dirJ.*potential;

            nDoF = obj.mdof.ndofs;
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
      obj.cbFreeze = false; % Unfreeze the branch of oedometric compressibility 

      % Get state and set the initial conditions
      state = getState(obj);
      state.strainAcc = state.strainAcc + state.strain;      
      mapSp = state.stress < state.Sp;
      state.Sp(mapSp) = state.stress(mapSp);
      state.e0(mapSp) = state.e(mapSp);

      % ------ Update the Sedimentation ------
      % Update the sedimentation 
      dt = state.time-getStateOld(obj,'time');
      sedAdd = dt*state.sedmRate;
      state.sedmAcc = state.sedmAcc + sedAdd;

      newCellSed = zeros(size(state.sedmAcc));
      mapNewCells(1:numel(obj.mdof.bathy),1) = false;
      
      % initialize the columns
      dzSedmAdd = sum(sedAdd,2);
      mapSedmPos = dzSedmAdd > 0;

      initCol = and(mapSedmPos,obj.initColumn);
      cellVar = and(mapSedmPos,~obj.initColumn);

      newCellSed(initCol,:) = sedAdd(initCol,:);
      mapNewCells(initCol) = true;
      
      if any(cellVar)
        % ------ Update the cells at the top (before the grow) ------
        dofs = getBordZMax(obj.mdof);
        idx = 1:length(dofs);
        dofs = dofs(cellVar);
        idx = idx(cellVar);

        dzOld = obj.cellDims(dofs,3);
        dzDef = dzOld.*(1+state.strain(cellVar));
        dzNew = dzDef+dzSedmAdd(cellVar);        

        growCols = dzNew > obj.heightControl*(1+obj.minCellHeight); % <-- Very important check
        
        state.cellDefm(dofs) = 0;
        state.strainAcc(dofs) = 0;
        if (any(~growCols))
          ngrowCols = ~growCols;
          dof = dofs(ngrowCols);
          loc = idx(ngrowCols);

          fracCell = dzOld(ngrowCols).*obj.matfrac(dof,:) + sedAdd(loc,:);
          obj.cellDims(dof,3) = dzNew(ngrowCols);
          obj.matfrac(dof,:) = fracCell./sum(fracCell,2);

          % Mean average
          % state.void(dof) = (dzDef(ngrowCols).*state.void(dof) +...
          %   dzSedmAdd(loc).*obj.voidTop(loc))./dzNew(ngrowCols);

          % Harmonic series
          state.void(dof) = ...
            (dzDef(ngrowCols) + dzSedmAdd(loc)) ./ ...
            ( dzDef(ngrowCols) ./ state.void(dof) + ...
            dzSedmAdd(loc) ./ obj.top_e(loc) );
          state.void0(dof) = ...
            (dzDef(ngrowCols) + dzSedmAdd(loc)) ./ ...
            ( dzDef(ngrowCols) ./ state.void0(dof) + ...
            dzSedmAdd(loc) ./ obj.top_e0(loc) );
        end

        if (any(growCols))
          dof = dofs(growCols);
          loc = idx(growCols);

          frac = (obj.heightControl-dzDef(growCols))./(dzSedmAdd(loc));

          fracCell = dzOld(growCols).*obj.matfrac(dof,:) + frac.*sedAdd(loc,:);
          obj.cellDims(dof,3) = obj.heightControl;
          obj.matfrac(dof,:) = fracCell./sum(fracCell,2);

          % Mean average
          % % % % state.void(dof) = (dzDef(growCols).*state.void(dof) +...
          % % % %   frac.*dzSedmAdd(loc).*obj.voidTop(loc))./obj.heightControl;

          % Harmonic series
          state.void(dof) = ...
            (dzDef(growCols) + frac.*dzSedmAdd(loc)) ./ ...
            ( dzDef(growCols) ./ state.void(dof) + ...
            (frac.*dzSedmAdd(loc)) ./ obj.top_e(loc) );

          state.void0(dof) = ...
            (dzDef(growCols) + frac.*dzSedmAdd(loc)) ./ ...
            ( dzDef(growCols) ./ state.void0(dof) + ...
            (frac.*dzSedmAdd(loc)) ./ obj.top_e0(loc) );

          newCellSed(loc,:) = (1-frac).*sedAdd(loc,:);
          mapNewCells(loc) = true;
        end
      end
      obj.initColumn(initCol) = false;
      
      state.strain(:) = 0;
      % ------ Grow the grid if the threshold is reached ------
      newcells = sum(mapNewCells);
      state.newcells = newcells;
      if newcells>0        
        ndofs = obj.mdof.ndofs;
        dofs = (ndofs+1:ndofs+newcells)';
        obj.matfrac(end+1:end+newcells,:) = newCellSed(mapNewCells,:)./sum(newCellSed(mapNewCells,:),2);

        % Update the grid and locate the new dofs
        newlayer = grow(obj.mdof,mapNewCells);
        % growLayer(obj.tmpObj,mapNewCells);
        updateConectivities(obj.tmpObj,mapNewCells,sum(newCellSed(mapNewCells,:),2));
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

        % Update the states
        state.maxDofUnchanged = getMaxDofUnchanged(obj.mdof);

        state.pressure(end+1:end+newcells) = 0.;
        state.stress(end+1:end+newcells) = 0.;
        state.strain(end+1:end+newcells) = 0.;
        state.strainAcc(end+1:end+newcells) = 0.;
        state.cellDefm(end+1:end+newcells) = 0.;
        state.Sp(end+1:end+newcells) = 0.;
        state.void(end+1:end+newcells) = 0.;
        state.void0(end+1:end+newcells) = 0.;
        state.age(end+1:end+newcells) = state.time;
        state.cb(end+1:end+newcells) = 0.;

        state.sedmAcc(mapNewCells,:) = newCellSed(mapNewCells,:);        
        state.void(dofs) = obj.top_e(mapNewCells);
        state.void0(dofs) = obj.top_e0(mapNewCells);

        frac = sum(newCellSed(mapNewCells,:),2)./sum(sedAdd(mapNewCells,:),2);
        stressAdd = - frac.*obj.top_dS(mapNewCells)*dt;

        stsInit = zeros(newcells,1);
        stsCellPreCons = zeros(newcells,1);
        matTmp = newCellSed(mapNewCells,:)./sum(newCellSed(mapNewCells,:),2);
        for mat=1:obj.nmat
          tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw;
          stsCellPreCons = stsCellPreCons + matTmp(:,mat).*tmpMat.Sp;
          stsInit = stsInit + matTmp(:,mat).*tmpMat.S0;
        end
        mapPreCon = stsCellPreCons<-stressAdd;
        stsCellPreCons(mapPreCon)=-stressAdd(mapPreCon);

        state.stress(dofs) = stressAdd;
        % state.stress(dofs) = stressAdd-stsInit;
        % state.Sp(dofs) = stressAdd;
        state.Sp(dofs) = -stsCellPreCons;
      end
      
      setStateOld(obj,state);
      setState(obj,state);
    end

    function updateState(obj,solution)
      % UPDATESTATE Iterate in the time step.
      state = getState(obj);
      stateOld = getStateOld(obj);

      % Some variables wrappers
      % Sprev = stateOld.stress;
      Scurr = state.stress + solution;
      Sp = stateOld.Sp;
      eprev = stateOld.e;

      % % New version
      map1 = Scurr<0;
      ep = state.e0;
      tmpMat = getCellsProp(obj,'sedmat');
      mapSS = SedimentMaterial.curveBranch(-Scurr(map1),-Sp(map1),tmpMat.S1(map1),tmpMat.S2(map1));
      ecurr = SedimentMaterial.computeVoid(-Scurr(map1),-Sp(map1),ep(map1),mapSS,tmpMat);

      % Update the Mesh Deformation - vertical deformation
      eps = (ecurr-eprev)./(1+eprev);
      % eps = delta_e./(1+eprev);

      % Update the Mesh Deformation - vertical deformation
      % eps = delta_e./(1+stateOld.void);

      % Lagrangian strain model
      % strain = stateOld.strain(dofs) + eps;  
      strain = eps;  
      dz = obj.mdof.dl(:,3)+state.cellDefm;

      % Eulerian strain model
      % strain = stateOld.strain(dofs) + eps./(1+eps);
      % dz = obj.cellDims(dofs,3) + state.cellDefm(dofs);

      % Update State
      state.stress = Scurr;
      state.pressure = state.pressure + solution;
      state.void = ecurr;
      % state.void = eprev+delta_e;
      state.strain = strain;
      state.cellDefm = stateOld.cellDefm + strain.*dz;

      obj.setState(state);
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
      p = getState(obj,obj.getField());
      pOld = getStateOld(obj,obj.getField());
      dt = getState(obj,'time') - getStateOld(obj,'time');

      rhsStiff = obj.H*p;
      rhsCap = (obj.P/dt)*(p - pOld);
      rhs = rhsStiff + rhsCap;

      % Computing the sedimentation contribution
      volCell = getCellVols(obj.mdof,getStateOld(obj,'cellDefm'));      
      cb = getState(obj,'cb');

      rhsSedm = volCell .* cb .* obj.mdof.top_dS(getMapFromDofs(obj.mdof))';
      rhs = rhs - rhsSedm;
    end

    function computeStiffMat(obj)
      % Compute the Stiffness Matrix
      ndofs = obj.mdof.ndofs;
      lw = 1/obj.domain.materials.getFluid().getSpecificWeight();

      [trans,neigh] = getFaceTransmissibility(obj.mdof,...
        getCellsProp(obj,'conductivity'),getStateOld(obj,'cellDefm'));
      trans= lw*trans;

      if ~isempty(trans)
        sumDiagTrans = accumarray(neigh(:),repmat(trans,[2,1]),[ndofs,1]);
        obj.H = sparse([neigh(:,1); neigh(:,2); (1:ndofs)'],...
          [neigh(:,2); neigh(:,1); (1:ndofs)'],...
          [-trans; -trans; sumDiagTrans], ndofs, ndofs);
      else
        obj.H = sparse(ndofs,ndofs);
      end
    end

    function computeCapMat(obj)
      % COMPUTECAPMAT Builds the storage matrix (P).
      %
      % Includes:
      %   - Rock compressibility
      %   - Fluid compressibility

      % Selecting the active dofs
      ndofs = obj.mdof.ndofs;
      dofs = 1:ndofs;

      % Computing the volumes
      volCell = getCellVols(obj.mdof,getStateOld(obj,'cellDefm'));
      % strain = getStateOld(obj,'strainAcc');
      % volCell = obj.mdof.vol.*(1-strain);

      % Computing the storage coefficient
      poro = getState(obj,'e');
      poro = poro./(1+poro);

      beta = obj.domain.materials.getFluid().getFluidCompressibility();
      oedoComp = getState(obj,'cb');
      PVal = (oedoComp+beta*poro).*volCell;

      % Create the P matrix.
      obj.P = sparse(dofs,dofs,PVal,ndofs,ndofs);
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      % Set mesh output.
      sed = sum(interpolate(obj.domain.state,fac,'sedmAcc'),2);

      obj.grid = makeMeshOutput(obj.mdof,obj.grid,obj.coordX,obj.coordY,...
        obj.coordZ,sed,obj.initColumn);
      outPrint = finalizeState(obj,fac);
      if isempty(obj.bath)
        outPrint.bath = zeros(size(outPrint.comp));
      else
        outPrint.bath = getBath(obj.mdof,obj.bath);
      end

      [cellData,pointData] = buildPrintStruct(obj,outPrint);
    end

    function writeSolution(obj,fac,tID)
      outPrint = finalizeState(obj,fac);
      obj.domain.outstate.results(tID).time = outPrint.time;
      obj.domain.outstate.results(tID).pressure = outPrint.pres;
      obj.domain.outstate.results(tID).porosity = outPrint.poro;
      obj.domain.outstate.results(tID).stress = outPrint.stress;
      obj.domain.outstate.results(tID).strain = outPrint.strain;
      obj.domain.outstate.results(tID).compaction = outPrint.comp;
    end

    function [cellStr,pointStr] = buildPrintStruct(obj,state)
      cellStr = repmat(struct('name', 1, 'data', 1), 6, 1);
      cellStr(1).name = 'overpressure';
      cellStr(1).data = state.pres;
      cellStr(2).name = 'stress';
      cellStr(2).data = state.stress;
      cellStr(3).name = 'strain';
      cellStr(3).data = state.strain;
      cellStr(4).name = 'void';
      cellStr(4).data = state.void;
      cellStr(5).name = 'porosity';
      cellStr(5).data = state.poro;
      cellStr(6).name = 'hydraulicConductivity';
      cellStr(6).data = state.cond;
      cellStr(7).name = 'age';
      cellStr(7).data = state.age;

      cellStr(8).name = 'oedometricCompresibility';
      cellStr(8).data = state.cb;

      cellStr(9).name = 'stressPreConsolidated';
      cellStr(9).data = state.Sp;

      celLast = length(cellStr);
      for mat=1:obj.nmat
        cellStr(mat+celLast).name = sprintf('material%02d', mat);
        cellStr(mat+celLast).data = obj.matfrac(:,mat);
      end

      if ~isempty(obj.bath)
        pointStr = repmat(struct('name', 1, 'data', 1), 2, 1);
      end
      pointStr(1).name = 'compaction';
      pointStr(1).data = state.comp;
      if ~isempty(obj.bath)
        pointStr(2).name = 'bathymetry';
        pointStr(2).data = state.bath;
      end

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
        if isfield(data.Grid,'xfile')
          obj.coordX = load(data.Grid.xfile);
        else
          obj.coordX = data.Grid.xcoord;
        end

        if isfield(data.Grid,'yfile')
          obj.coordY = load(data.Grid.yfile);
        else
          obj.coordY = data.Grid.ycoord;
        end

        if isfield(data.Grid,'zfile')
          obj.coordZ = load(data.Grid.zfile);
        else
          obj.coordZ = data.Grid.zcoord;
        end
      else
        gresLog().error("The grid parameters for the simulation is not well defined!");
      end
      obj.heightControl = data.NewCellHeightControl;
    end

    function prepareMaterialFractions(obj,data)
      % Build the material fractions for each cell
      ndofs = obj.mdof.ndofs;
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
    end

    function computeCb(obj)
      % Some variables
      state = getState(obj);
      % stateOld = getStateOld(obj);

      % Some variables wrappers - New Version
      Scurr= state.stress;
      % Sp   = stateOld.Sp;
      void = state.e;

      map1 = Scurr<0;
      oedoComp = zeros(length(void),1);
      tmpMat = getCellsProp(obj,'sedmat');
      oedoComp(map1) = SedimentMaterial.computeOedo(-Scurr(map1),...
        void(map1),obj.cellMatBranch(map1,:),tmpMat);
      setState(obj,oedoComp,'cb');
    end

    function MaterialMap(obj)
      if ~obj.cbFreeze
        Scurr = getState(obj,'stress');
        Sp = getState(obj,'Sp');
        obj.cellMatBranch = false(obj.mdof.ndofs,4);
        map1 = Scurr<0;
        tmpMat = getCellsProp(obj,'sedmat');
        obj.cellMatBranch(map1,:) = SedimentMaterial.curveBranch(-Scurr(map1),-Sp(map1),tmpMat.S1(map1),tmpMat.S2(map1));
      end
    end

    function out = getCellsProp(obj,type,dofs)
      if ~exist("dofs","var")
        dofs = 1:obj.mdof.ndofs;
      end
      matfrac = obj.mdof.matfrac(dofs,:);

      type = lower(type);
      switch type
        case 'conductivity'
          out = zeros(length(dofs),3);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getConductivity();
            out = out + matfrac(:,mat).*tmpMat;
          end
          % try to impose kz0*10^((e-e0)/Ck)
          % dvoid = obj.getState("void")-obj.void0;
          % out(:,3)=out(:,3).*10.^(dvoid/10);
          % out(:,1:2)=[3*out(:,3) 3*out(:,3)];
        case 'compressidx'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getCompressibilityIdx();
            out = out + matfrac(:,mat).*tmpMat;
          end
        case 'recompressidx'
          out = zeros(length(dofs),1);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getReCompressibilityIdx();
            out = out + matfrac(:,mat).*tmpMat;
          end
        case 'sedmat'
          out=struct();
          out.Cc = zeros(obj.mdof.ndofs,1);
          out.Cr = zeros(obj.mdof.ndofs,1);
          out.S1 = zeros(obj.mdof.ndofs,1);
          out.S2 = zeros(obj.mdof.ndofs,1);
          out.emin = zeros(obj.mdof.ndofs,1);
          out.e1 = zeros(obj.mdof.ndofs,1);
          out.e2 = zeros(obj.mdof.ndofs,1);
          out.cbmin = zeros(obj.mdof.ndofs,1);
          out.kDecay = zeros(obj.mdof.ndofs,1);

          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw;
            out.Cc = out.Cc + matfrac(:,mat).*tmpMat.Cc;
            out.Cr = out.Cr + matfrac(:,mat).*tmpMat.Cr;
            out.S1 = out.S1 + matfrac(:,mat).*tmpMat.S1;
            out.S2 = out.S2 + matfrac(:,mat).*tmpMat.S2;
            out.emin = out.emin + matfrac(:,mat).*tmpMat.emin;
            out.e1 = out.e1 + matfrac(:,mat).*tmpMat.e1;
            out.e2 = out.e2 + matfrac(:,mat).*tmpMat.e2;
            out.cbmin = out.cbmin + matfrac(:,mat).*tmpMat.cbmin;
            out.kDecay = out.kDecay + matfrac(:,mat).*tmpMat.kDecay;
          end
        % case 'voidlowerlimit'
        %   out = zeros(length(dofs),1);
        %   for mat=1:obj.nmat
        %     tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw.getVoidLowerLim();
        %     out = out + obj.matfrac(dofs,mat).*tmpMat;
        %   end
        otherwise
          out = [];
      end
    end

    function data = finalizeState(obj,fac)
      % append state variable to output structure
      state = interpolate(obj.domain.state,fac);

      data.time = state.time;
      data.pres = state.pressure;
      data.stress = -state.stress;
      data.Sp = -state.Sp;
      data.strain = -state.strainAcc;
      data.void = state.void;
      data.age = getState(obj,'time')-state.age;
      data.cond = getCellsProp(obj,'conductivity');
      data.poro = state.void./(1+state.void);
      data.cb = state.cb;
      data.comp = -getComp(obj.mdof,state.cellDefm);
    end

  end

  methods (Static)

    function [void, stress] = initialCellProp(voidInit,stsParcial,stsInit,Cc,tolerance,maxNIter)
      % INITIALCELLPROP function to initialize the cell void rate and stress
      % using a iterative procedure.
      %
      % $$ \sigma_{trial} = (\frac{1}{1+e})*\sigma_{parcial} $$,
      %
      % $$ e_{c} = e_{0} - C_{r} \log(\frac{\sigma_{trial}}{\sigma_{init}}) $$
      %
      % where:
      % - $$e_{0}$$ - the initial void rate
      % - $$e_{c}$$ - the current void rate
      % - $$C_{r}$$ - the compression index
      % - $$\sigma_{init}$$ - the initial stress
      % - $$\sigma_{trial}$$ - the trial stress
      %
      % The iterative process consist in evaluate the stress using the void
      % rate and the void rate using the stress. This procedure is maintain
      % util the difference of void rate is smaller than a tolerance value
      % between iterations.
            
      % First iteration
      stress = (1./(1+voidInit)).*stsParcial;
      void = SedimentMaterial.getVoidRatioFromRef(stress,stsInit,voidInit,Cc);
      error = max(abs((void-voidInit)./voidInit));
      voidIt1 = void;

      % Iterative Procedure
      iter = 2;
      while and(error > tolerance, iter<maxNIter)
        stress = (1./(1+void)).*stsParcial;
        void = SedimentMaterial.getVoidRatioFromRef(stress,stsInit,voidInit,Cc);
        error = max(abs((void-voidInit)./voidInit));
        iter = iter+1;
      end

      if (iter>maxNIter)
        gresLog().warning(-1,'The iterative procedure for the void rate did not converge!');
      end

      hasIm = imag(void) ~= 0;
      if any(hasIm)
        void(hasIm) = voidIt1(hasIm);
        void = real(void);
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