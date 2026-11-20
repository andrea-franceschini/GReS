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

    matfrac (:,:)                % Material fractions per cell

    cellMatBranch
    cbFreeze = false

    deltaStress (:,1)
    voidTop
    e0Top

    initColumn (:,1)  % indication of the necessity of initialize the column
    bath (:,:)
  end

  properties (Access = protected)
    fieldId
  end

  properties (Access = private)
    nmat
    tol = 1e-8;
    niterMx = 100;
    minCellHeight = 1e-9;
    void0MixtureType = true;
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
      ncols = prod(map(1:2));
      segmX = diff(obj.coordX);
      segmY = diff(obj.coordY);
      segmZ = diff(obj.coordZ);
      [idI,idJ,idK] = getIJKFromDofs(obj.mdof);
      obj.cellDims = [segmX(idI) segmY(idJ) segmZ(idK)];      
      obj.initColumn(1:ncols,1) = true;

      % open space for sediments
      obj.deltaStress = zeros(ncols,1);
      obj.voidTop = zeros(ncols,1);
      obj.e0Top = zeros(ncols,1);

      % Build the material fractions for each cell
      prepareMaterialFractions(obj,input.Domain.Initial);

      % Prepare the faces connectivity
      nelm = obj.mdof.ndofs;
      dofs = (1:nelm)';      
      cellNeigh = getNeigh(obj.mdof);
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

      state.pressure = zeros(nelm,1);
      state.stress = zeros(nelm,1);
      state.strain = zeros(nelm,1);
      state.strainAcc = zeros(nelm,1);
      state.cellDefm = zeros(nelm,1);
      state.Sp = zeros(nelm,1);      
      state.void = zeros(nelm,1);
      state.age = zeros(nelm,1);
      state.void0 = zeros(nelm,1);
      state.cb = zeros(nelm,1);

      state.sedmRate = zeros(prod(map(1:2)),obj.nmat);
      state.sedmAcc  = zeros(prod(map(1:2)),obj.nmat);

      setState(obj,state)

      % Set bathymetry
      if isfield(input,"Bathymetry")
        vls = readmatrix(input.Bathymetry.file,'Delimiter','\t');
        mapDim = size(vls);
        if all(mapDim==map(1:2))
          tmp = vls;
        else
          [x_ref, y_ref] = meshgrid(0:1/(mapDim(1)-1):1, 0:1/(mapDim(2)-1):1);
          [x_new, y_new] = meshgrid(0:1/(map(1)-1):1, 0:1/(map(2)-1):1);
          tmp = interp2(x_ref, y_ref, vls, x_new, y_new, 'cubic')';
        end
        obj.bath = tmp;
      end
    end

    function initialize(obj)
      % Initialize the Sedimentation solver before the simulation starts

      % Initialize the States
      state = getState(obj);
      state.age(:) = state.time;
      ncells = [length(obj.coordX)-1 length(obj.coordY)-1 length(obj.coordZ)-1];

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
          stsCellPreCons = zeros(ncols,1);
          Cr = zeros(ncols,1);
          Cc = zeros(ncols,1);
          for mat=1:obj.nmat
            tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw;

            DStress = DStress + (tmpMat.gamma - gamma_w)*matFract(:,mat).*dh;
            voidTmp = voidTmp + matFract(:,mat).*tmpMat.e0;
            Cc = Cc + matFract(:,mat).*tmpMat.Cc;
            Cr = Cr + matFract(:,mat).*tmpMat.Cr;
            stsCellAcc = stsCellAcc + matFract(:,mat).*tmpMat.S0;
            stsCellPreCons = stsCellPreCons + matFract(:,mat).*tmpMat.Sp;
          end

          % Iterative procedure
          strial = (1./(1+voidTmp)).*DStress;
          % strial = (1./(1+voidTmp)).*DStress+stsCons(map);
          [void,~] = Sedimentation.initialCellProp(voidTmp,-strial,...
            -stsCellAcc,-stsCons(map),Cr,obj.tol,obj.niterMx);

          stress = (1./(1+void)).*DStress;
          % stress = (1./(1+void)).*DStress+stsCons(map);

          mapPreCon = stsCellPreCons<stress;
          stsCellPreCons(mapPreCon)=stress(mapPreCon);

          stStore(map) = stress;
          state.stress(dof) = -stress;
          % state.stress(dof) = -stress-1;          
          % state.Sp(dof) = -stress;
          state.Sp(dof) = -stsCellPreCons;
          state.void(dof) = void;

          state.void0 = SedimentMaterial.getVoidPreCon(stress,stsCellPreCons,void,Cr,Cc);
        end
      end

      % Set States
      obj.setState(state);
      obj.setStateInit(state);
      obj.setStateOld(state);
    end

    function timeStepSetup(obj)
      % initialize the physics solver for the time step
      dofs = getDofs(obj.mdof);
      state = getState(obj);
      stateOld = getStateOld(obj);

      ncols = prod(obj.mdof.ncells(1:2));
      obj.deltaStress(:) = 0.;
      obj.voidTop(:) = 0.;
      obj.e0Top(:) = 0.;

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
        obj.voidTop = obj.voidTop + matFract(:,mat).*tmpMat.e0;
        Cc = Cc + matFract(:,mat).*tmpMat.Cc;
        Cr = Cr + matFract(:,mat).*tmpMat.Cr;
        stsInit = stsInit + matFract(:,mat).*tmpMat.S0;
      end
      % obj.deltaStress = sedmWeight;

      % iterate to find the true value of e and stress
      % stsCons = zeros(ncols,1);
      strial = (1./(1+obj.voidTop(map))).*sedmWeight(map);
      [void,~] = Sedimentation.initialCellProp(obj.voidTop(map),-strial,...
        -stsInit(map),-stsInit(map),Cr(map),obj.tol,obj.niterMx);

      obj.deltaStress(map) = (1./(1+void)).*sedmWeight(map);
      obj.voidTop(map) = void;
      obj.e0Top(map) = SedimentMaterial.getVoidPreCon(obj.deltaStress(map),stsInit(map),void,Cr(map),Cc(map));

      dstress = obj.deltaStress(getMapFromDofs(obj.mdof,dofs));
      state.stress(dofs) = state.stress(dofs) - dstress*dt;

      setState(obj,state);

      % Initialize the cells cb branches
      obj.cellMatBranch = false(numel(state.stress),4);
      map = state.stress < 0;
      tmpMat = getCellsProp(obj,'sedmat');
      obj.cellMatBranch(map,:) = SedimentMaterial.curveBranch(...
        -state.stress(map),-state.Sp(map),tmpMat.S1(map),tmpMat.S2(map));

      if obj.cbFreeze
        gresLog().log(0,' %s \n','-----------------------------------------------------------');
        gresLog().log(0,' %s \n','IN THE NEXT TIME STEP, THE CB CURVE IS FREEZED!');
      end
    end

    function assembleSystem(obj,dt)
      % Update the cells half transmissibility
      obj.MaterialMap();
      obj.computeHalfTrans();
      obj.computeCb();

      % Assembling the system.
      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);
    end

    function applyBC(obj,bcId,t)
      for bcId = obj.bcs
        dofs = obj.mdof.getBord(bcId.surface);
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
            trans = obj.halfTrans(dofs,axis);
            mu = obj.domain.materials.getFluid().getSpecificWeight();
            dirJ = 1/mu*trans;
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
      % ADVANCESTATE Finalizes the time step and updates the grid.
      obj.cbFreeze = false;
      state = getState(obj);

      % Update accumulated strain and pre consolidation state
      state.strainAcc = state.strainAcc + state.strain;
      mapSp = state.stress < state.Sp;
      state.Sp(mapSp) = state.stress(mapSp);
      state.void0(mapSp) = state.void(mapSp);

      % Sediment deposited during the current time step
      dt = state.time - getStateOld(obj,'time');
      sedAdd = dt*state.sedmRate;
      state.sedmAcc = state.sedmAcc + sedAdd;

      % Update the growing cells at the top of the domain
      [state,sedNewCells] = updateGrowingCells(obj,state,sedAdd);

      % Reset the strain deformation in all domain for the next time step.
      state.strain(:) = 0;

      % Create and initialize new cells
      state.newcells = nnz(any(sedNewCells,2));
      if state.newcells > 0
        [state,newDofs] = growGrid(obj,state,sedNewCells);
        state = initializeStatesNewCells(obj,state,newDofs,sedNewCells,sedAdd);
      end

      % Commit the updated state
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
      eprev = stateOld.void;

      % % New version
      map1 = Scurr<0;
      ep = state.void0;
      tmpMat = getCellsProp(obj,'sedmat');
      % ecurr = SedimentMaterial.computeVoid(-Scurr(map1),-Sp(map1),ep(map1),obj.cellMatBranch(map1,:),tmpMat);
      mapSS = SedimentMaterial.curveBranch(-Scurr(map1),-Sp(map1),tmpMat.S1(map1),tmpMat.S2(map1));
      ecurr = SedimentMaterial.computeVoid(-Scurr(map1),-Sp(map1),ep(map1),mapSS,tmpMat);

      % % % % Old version
      % % % Cc = getCellsProp(obj,'compressIdx');
      % % % Cr = getCellsProp(obj,'recompressIdx');
      % % % delta_e = SedimentMaterial.getDeltaVoidRatio(Scurr,Sprev,Sp,Cc,Cr);
      % % % ecurr = eprev + delta_e;

      % Update the Mesh Deformation - vertical deformation
      eps = (ecurr-eprev)./(1+eprev); % - New version
      % eps = delta_e./(1+eprev);

      % Update the Mesh Deformation - vertical deformation
      % eps = delta_e./(1+stateOld.void);

      % Lagrangian strain model
      % strain = stateOld.strain(dofs) + eps;  
      strain = eps;  
      dz = obj.cellDims(:,3)+state.cellDefm;

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
      
      dl = getStateOld(obj,'cellDefm');
      volCell = prod(obj.cellDims(:,1:2),2).*(obj.cellDims(:,3) + dl);
      cb = getState(obj,'cb');

      rhsSedm = volCell .* cb .* obj.deltaStress(getMapFromDofs(obj.mdof));
      rhs = rhs - rhsSedm;
    end

    function computeStiffMat(obj)
      % Compute the Stiffness Matrix
      ndofs = obj.mdof.ndofs;
      lw = 1/obj.domain.materials.getFluid().getSpecificWeight();

      idi = sub2ind([ndofs 3], obj.facesNeigh(:,1), obj.facesNeighDir(:));
      idk = sub2ind([ndofs 3], obj.facesNeigh(:,2), obj.facesNeighDir(:));
      Tii = obj.halfTrans(idi);
      Tik = obj.halfTrans(idk);
      trans = lw./(1./Tii+1./Tik);

      sumDiagTrans = accumarray(obj.facesNeigh(:),repmat(trans,[2,1]),[ndofs,1]);
      obj.H = sparse([obj.facesNeigh(:,1); obj.facesNeigh(:,2); (1:ndofs)'],...
        [obj.facesNeigh(:,2); obj.facesNeigh(:,1); (1:ndofs)'],...
        [-trans; -trans; sumDiagTrans], ndofs, ndofs);
    end

    function computeCapMat(obj)
      % COMPUTECAPMAT Builds the storage matrix (P).
      %
      % Includes:
      %   - Rock compressibility
      %   - Fluid compressibility

      % Selecting the active dofs
      ndofs = obj.mdof.ndofs;
      dofs = getDofs(obj.mdof);

      % Computing the volumes
      dl = getStateOld(obj,'cellDefm');
      volCell = prod(obj.cellDims(:,1:2),2).*(obj.cellDims(:,3) + dl);

      % Computing the storage coefficient
      poro = getState(obj,'void');
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

    function computeHalfTrans(obj)
      % COMPUTEHALFTRANS Computes cells half-transmissibility.
      obj.halfTrans = zeros(obj.mdof.ndofs,3);

      % Get mesh dimension and update with the deformation
      dl = getStateOld(obj,'cellDefm');
      dx = obj.cellDims(:,1);
      dy = obj.cellDims(:,2);
      dz = obj.cellDims(:,3) + dl;

      % Computing half transmissibilities
      condCell = getCellsProp(obj,'conductivity');
      obj.halfTrans(:,1)= (dy.*dz)./(dx/2).*condCell(:,1);
      obj.halfTrans(:,2)= (dx.*dz)./(dy/2).*condCell(:,2);
      obj.halfTrans(:,3)= (dx.*dy)./(dz/2).*condCell(:,3);
    end

    function computeCb(obj)
      % Some variables
      state = getState(obj);
      % stateOld = getStateOld(obj);

      % Some variables wrappers - New Version
      Scurr= state.stress;
      % Sp   = stateOld.Sp;
      void = state.void;

      map1 = Scurr<0;
      oedoComp = zeros(length(void),1);
      tmpMat = getCellsProp(obj,'sedmat');
      oedoComp(map1) = SedimentMaterial.computeOedo(-Scurr(map1),void(map1),obj.cellMatBranch(map1,:),tmpMat);

      % % Some variables wrappers - Old Version
      % Sprev = stateOld.stress;
      % Scurr = state.stress;
      % Sp = stateOld.Sp;
      % void = state.void;
      % Cc = getCellsProp(obj,'compressIdx');
      % Cr = getCellsProp(obj,'recompressIdx');
      % oedoComp = (1./(1+void)).*SedimentMaterial.getDevVoidRatio(Scurr,Sprev,Sp,Cc,Cr);

      setState(obj,oedoComp,'cb');
    end

    function MaterialMap(obj)
      if ~obj.cbFreeze
        Scurr = getState(obj,'stress');
        Sp = getState(obj,'Sp');
        obj.cellMatBranch = false(length(Scurr),4);
        map1 = Scurr<0;
        tmpMat = getCellsProp(obj,'sedmat');
        obj.cellMatBranch(map1,:) = SedimentMaterial.curveBranch(-Scurr(map1),-Sp(map1),tmpMat.S1(map1),tmpMat.S2(map1));
      end
    end


    % function for deal with the advance state
    function [state,sedNewCells] = updateGrowingCells(obj,state,sedAdd)

      sedNewCells = zeros(size(sedAdd));
      dzSed = sum(sedAdd,2);
      mapSed = dzSed > 0;

      % Initialize empty columns
      initCol = mapSed & obj.initColumn;
      sedNewCells(initCol,:) = sedAdd(initCol,:);
      obj.initColumn(initCol) = false;

      % Existing active columns receiving sediment
      cellVar = mapSed & ~obj.initColumn & ~initCol;

      if ~any(cellVar)
        return
      end

      dofs = getBordZMax(obj.mdof);
      dofs = dofs(cellVar);
      colIds = find(cellVar);

      % Current deformed geometry
      % dzOld = obj.cellDims(dofs,3);
      dzDef = obj.cellDims(dofs,3) + state.cellDefm(dofs);
      dzAdd = dzSed(colIds);

      dzTotal = dzDef + dzAdd;
      Hc = obj.heightControl;

      % Determine which columns require division
      growCols = dzTotal > Hc*(1+obj.minCellHeight);

      % Existing cell height after deposition
      dzLower = dzTotal;

      % For growing columns, preserve the existing deformed height if it
      % already exceeds Hc.
      dzLower(growCols) = max(Hc,dzDef(growCols));

      % Sediment incorporated into existing cells
      dzIncorporated = max(0,dzLower - dzDef);

      % Fraction of deposited sediment incorporated
      frac = dzIncorporated./dzAdd;

      % Update existing cell properties      
      state = updateGrowingCellsProperties(obj,state,dofs,colIds,dzDef, ...
        dzLower,dzIncorporated,sedAdd);

      % Update the reference geometry
      obj.cellDims(dofs,3) = dzLower;

      % Deformation is now included in cellDims
      state.cellDefm(dofs) = 0;
      state.strainAcc(dofs) = 0;

      % Sediment assigned to new cells
      if any(growCols)
        loc = colIds(growCols);
        sedNewCells(loc,:) = (1-frac(growCols)).*sedAdd(loc,:);
      end
    end

    function state = updateGrowingCellsProperties(obj,state,dofs,colIds,dzDef,dzLower,dzInc,sedAdd)

      % Avoid mixing when no sediment is incorporated
      active = dzInc > 0;
      if ~any(active)
        return
      end

      % Original cell properties
      eOld = state.void(dofs);
      matOld = obj.matfrac(dofs,:);

      % Newly deposited sediment properties
      eSed = obj.voidTop(colIds);

      dof = dofs(active);
      loc = colIds(active);

      hOld = dzDef(active);
      hSed = dzInc(active);
      hTotal = dzLower(active);

      e1 = eOld(active);
      e2 = eSed(active);

      % Solid-volume-conserving void ratio
      solidOld = hOld./(1+e1);
      solidSed = hSed./(1+e2);
      state.void(dof) = hTotal./(solidOld+solidSed)-1;

      % Sediment material fractions
      sedFrac = sedAdd(loc,:)./sum(sedAdd(loc,:),2);

      % Solid-volume-weighted material fractions
      matMix = solidOld.*matOld(active,:) + solidSed.*sedFrac;
      obj.matfrac(dof,:) = matMix./sum(matMix,2);

      % How to recompute the pre consolidated void rate
      if obj.void0MixtureType
        % Solid-volume-conserving void ratio for the pre consolidated level

        % Original cell properties
        eOld = state.void0(dofs);
        eSed = obj.e0Top(colIds);
        e1 = eOld(active);
        e2 = eSed(active);

        % Solid-volume-conserving void ratio
        solidOld = hOld./(1+e1);
        solidSed = hSed./(1+e2);
        state.void0(dof) = hTotal./(solidOld+solidSed)-1;
      else
        % Recompute the pre consolidation void ratio for the new cell mixture
        voidMix = state.void(dof);
        fracMix = obj.matfrac(dof,:);

        % Initialize equivalent material properties and effective stress
        CcMix = zeros(size(voidMix));
        CrMix = zeros(size(voidMix));
        SpMix = zeros(size(voidMix));
        stressMix = zeros(size(voidMix));
        gammaW = obj.domain.materials.getFluid().getSpecificWeight();

        for mat = 1:obj.nmat
          matLaw = obj.domain.materials.getMaterial(mat).ConstLaw;
          CcMix = CcMix + fracMix(:,mat).*matLaw.Cc;
          CrMix = CrMix + fracMix(:,mat).*matLaw.Cr;
          SpMix = SpMix + fracMix(:,mat).*matLaw.Sp;
          stressMix = stressMix + (matLaw.gamma - gammaW).*fracMix(:,mat).*hTotal;
        end
        stressMix = stressMix./(1+voidMix);

        % Update preconsolidation void ratio
        state.void0(dof) = SedimentMaterial.getVoidPreCon(stressMix, SpMix, voidMix, CrMix, CcMix);
      end
    end

    function [state,dofs] = growGrid(obj,state,sediments)
      map = any(sediments,2);
      ncells = state.newcells;
      ndofs = obj.mdof.ndofs;
      dofs = (ndofs+1:ndofs+ncells)';

      % New cell material fractions
      sedNew = sediments(map,:);
      dzNew = sum(sedNew,2);
      obj.matfrac(end+1:end+ncells,:) = sedNew./dzNew;

      % Grow the DOF manager
      newlayer = grow(obj.mdof,map);

      if newlayer
        obj.coordZ(end+1) = obj.coordZ(end)+obj.heightControl;
      end

      % New cell dimensions
      segmX = diff(obj.coordX);
      segmY = diff(obj.coordY);
      dx = repmat(segmX,[length(segmY) 1]);
      dy = repelem(segmY,length(segmX));
      obj.cellDims(end+1:end+ncells,:) = [dx(map),dy(map),dzNew];

      % Update face connectivity
      updateConnectivity(obj,dofs,map);
    end

    function updateConnectivity(obj,dofs,map)
      newcells = numel(dofs);

      nlaysByCol = obj.mdof.laysByCol(:);
      atTop = nlaysByCol == obj.mdof.ncells(3);

      newCellNotTop = ~atTop(map);

      neigh = getNeigh(obj.mdof,dofs);

      face = zeros(5*newcells,2);
      faceAxis = ones(5*newcells,1);
      faceAct = false(5*newcells,1);

      % Vertical connections
      face(1:newcells,:) = [neigh(:,5),dofs];
      faceAxis(1:newcells) = 3;
      faceAct(1:newcells) = true;

      % Horizontal connections
      for ref = 1:4
        tmp = neigh(:,ref);
        flag0 = (tmp < dofs) & (tmp ~= 0) & newCellNotTop;
        flag1 = tmp > dofs;
        flag = flag0 | flag1;
        nfaces = nnz(flag);

        interval = ref*newcells + (1:nfaces);

        face(interval,:) = [dofs(flag),tmp(flag)];
        faceAxis(interval) = ceil(ref/2);
        faceAct(interval) = true;
      end

      % Append new faces
      obj.facesNeigh = [ obj.facesNeigh; face(faceAct,:)];
      obj.facesNeighDir = [obj.facesNeighDir; faceAxis(faceAct) ];
    end

    function state = initializeStatesNewCells(obj,state,dofs,newCellSed,sedAdd)
      newcells = state.newcells;

      % Extend state variables
      state.pressure(end+1:end+newcells) = 0;
      state.stress(end+1:end+newcells) = 0;
      state.strain(end+1:end+newcells) = 0;
      state.strainAcc(end+1:end+newcells) = 0;
      state.cellDefm(end+1:end+newcells) = 0;
      state.Sp(end+1:end+newcells) = 0;
      state.void(end+1:end+newcells) = 0;
      state.void0(end+1:end+newcells) = 0;
      state.age(end+1:end+newcells) = state.time;
      state.cb(end+1:end+newcells) = 0;

      map = any(newCellSed,2);

      % Newly deposited sediment
      sedNew = newCellSed(map,:);

      % Preserve the original top-cell accumulation convention
      state.sedmAcc(map,:) = sedNew;

      % Initial void properties
      state.void(dofs) = obj.voidTop(map);
      state.void0(dofs) = obj.e0Top(map);

      % Fraction of the current deposition assigned to new cells
      frac = sum(sedNew,2)./sum(sedAdd(map,:),2);

      % Initial effective stress
      StressNewCells = -frac.*obj.deltaStress(map);

      % Initial preconsolidation stress
      SpNewCells = zeros(newcells,1);
      matFrac = sedNew./sum(sedNew,2);
      for mat = 1:obj.nmat
        tmpMat = obj.domain.materials.getMaterial(mat).ConstLaw;
        SpNewCells = SpNewCells - matFrac(:,mat).*tmpMat.Sp;
      end
      SpNewCells = min(SpNewCells,StressNewCells);

      state.stress(dofs) = StressNewCells;
      state.Sp(dofs) = SpNewCells;
    end



    function out = getCellsProp(obj,type,dofs)
      if ~exist("dofs","var")
        dofs = getDofs(obj.mdof);
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
          % dvoid = obj.getState("void")-obj.void0;
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
            out.Cc = out.Cc + obj.matfrac(dofs,mat).*tmpMat.Cc;
            out.Cr = out.Cr + obj.matfrac(dofs,mat).*tmpMat.Cr;
            out.S1 = out.S1 + obj.matfrac(dofs,mat).*tmpMat.S1;
            out.S2 = out.S2 + obj.matfrac(dofs,mat).*tmpMat.S2;
            out.emin = out.emin + obj.matfrac(dofs,mat).*tmpMat.emin;
            out.e1 = out.e1 + obj.matfrac(dofs,mat).*tmpMat.e1;
            out.e2 = out.e2 + obj.matfrac(dofs,mat).*tmpMat.e2;
            out.cbmin = out.cbmin + obj.matfrac(dofs,mat).*tmpMat.cbmin;
            out.kDecay = out.kDecay + obj.matfrac(dofs,mat).*tmpMat.kDecay;
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

    function [void, stress] = initialCellProp(voidInit,stsParcial,stsInit,stsCons,Cc,tolerance,maxNIter)
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
      % - $$C_{r}$$ - the compression index
      % - $$\sigma_{init}$$ - the initial stress
      % - $$\sigma_{trial}$$ - the trial stress
      % - $$cte$$ - a constant stress value
      %
      % The iterative process consist in evaluate the stress using the void
      % rate and the void rate using the stress. This procedure is maintain
      % util the difference of void rate is smaller than a tolerance value
      % between iterations.
            
      % First iteration
      stress = (1./(1+voidInit)).*stsParcial+stsCons;
      void = SedimentMaterial.getVoidRatioFromRef(stress,stsInit,voidInit,Cc);
      error = max(abs((void-voidInit)./voidInit));
      voidIt1 = void;

      % Iterative Procedure
      iter = 2;
      while and(error > tolerance, iter<maxNIter)
        stress = (1./(1+void)).*stsParcial+stsCons;
        void = SedimentMaterial.getVoidRatioFromRef(stress,stsInit,voidInit,Cc);
        error = max(abs((void-voidInit)./voidInit));
        iter=iter+1;
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