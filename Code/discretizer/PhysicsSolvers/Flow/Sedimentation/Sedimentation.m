classdef Sedimentation < PhysicsSolver
  % Sedimentation model subclass

  properties
     mapSediments = SedimentsMap()

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

      if ~(isfield(input,'domain') || isfield(input,'Domain'))
        error("Domain for the simulation not defined!");
      end
      obj.nmat = length(obj.materials.db)-1;
      % obj.mesh = gridForSedimentation( "XML",input.domain, ...
      %   "Materiais",obj.materials.matMap);
      obj.mesh = gridForSedimentation( "XML",input.domain, ...
        "NumMateriais",obj.nmat);

            
      % Initialize the sedimentation control
      obj.mapSediments = SedimentsMap(obj.mesh.ncells(1:2),input.sediment_map);
      % obj.mapSediments = SedimentsMap(obj.grid.ncells(1:2),input.sediment_map);

      % Initialize the BCs
      obj.prepareBC(input.boundary);

      % Initialize the states
      obj.getState().data.(obj.getField()) = zeros(obj.mesh.getNumberCells,1);

      % Initialize the Transmissibility.
      prepareSystem(obj);

      obj.fieldId = 1;

      obj.domain.J{obj.fieldId,obj.fieldId} = [];
      obj.domain.rhs{obj.fieldId} = [];

      % Compute the initial system.
      computeHPInitial(obj);


    end

    function assembleSystem(obj,dt)
      obj.domain.J{obj.fieldId,obj.fieldId} = computeMat(obj,dt);
      obj.domain.rhs{obj.fieldId} = computeRhs(obj,dt);
    end

    function applyBC(obj,bcId,t) % <---------HERE
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

    function updateState(obj,solution)
      obj.Solver.updateState(solution);
    end

    function [cellData,pointData] = writeVTK(obj,fac,t)
      [cellData,pointData] = obj.Solver.writeVTK(fac,t);
    end

    function writeMatFile(obj,fac,tID)
      obj.Solver.writeMatFile(fac,tID);
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
      [cellId, faceArea, dz] = obj.mesh.getBordCell(bc.surface);
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
        case 'newmann'
          v=bc.data.*ones(length(cellId),1);
          vals = sum(vecN.*v,2);
        case 'dirichlet'
          gamma = obj.materials.getFluid().getFluidSpecWeight();
          mu = obj.materials.getFluid().getDynViscosity();
          permCell=0;
          for mat=1:obj.nmat
            tmpMat=obj.materials.getMaterial(mat).PorousRock.getPermVoigt();
            permCell = permCell+obj.mesh.matfrac(cellId,mat)*tmpMat(axis);
          end
          dirJ = 1/mu*(faceArea.*permCell);
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
        %   dz = obj.mesh.cellCentroid(cellId,3) - obj.faces.faceCentroid(faceID,3);
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
      ents = obj.mesh.mapCellIds(:);
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
      ncells = obj.mesh.numberActiveCells;
      tmpVec = lw.*obj.trans;
      sumDiagTrans = accumarray(obj.facesNeigh(:),repmat(tmpVec,[2,1]),[ncells,1]);
      obj.H = sparse([obj.facesNeigh(:,1); obj.facesNeigh(:,2); (1:ncells)'],...
        [obj.facesNeigh(:,2); obj.facesNeigh(:,1); (1:ncells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], ncells, ncells);
    end

    function computeInitialCapMat(obj)
      ncells = obj.mesh.numberActiveCells;
      beta = obj.materials.getFluid().getFluidCompressibility();
      PVal = obj.rockCmCell+beta*obj.poroCell;
      obj.P = PVal.*speye(ncells);

      % subCells = obj.dofm.getFieldCells(obj.fieldId);
      % nSubCells = length(subCells);
      % poroMat = zeros(nSubCells,1);
      % alphaMat = zeros(nSubCells,1);
      % beta = obj.materials.getFluid().getFluidCompressibility();
      % for m = 1:obj.mesh.nCellTag
      %   if ~ismember(m,obj.dofm.getTargetRegions([obj.getField(),"displacements"]))
      %     % compute alpha only if there's no coupling in the
      %     % subdomain
      %     alphaMat(m) = obj.materials.getMaterial(m).ConstLaw.getRockCompressibility();
      %   end
      %   poroMat(m) = obj.materials.getMaterial(m).PorousRock.getPorosity();
      % end
      % % (alpha+poro*beta)
      % PVal = alphaMat(obj.mesh.cellTag(subCells)) + beta*poroMat(obj.mesh.cellTag(subCells));
      % PVal = PVal.*obj.mesh.cellVolume(subCells);
      % nDoF = obj.dofm.getNumbDoF(obj.fieldId);
      % [~,~,dof] = unique(subCells);
      % obj.P = sparse(dof,dof,PVal,nDoF,nDoF);
    end

    function gTerm = finalizeRHSGravTerm(obj,lw)
      gTerm = accumarray(obj.facesNeigh(:), ...
        [lw.*obj.rhsGrav; -lw.*obj.rhsGrav],[obj.mesh.numberActiveCells,1]);
      gTerm = gTerm(obj.mesh.getActiveCells);
    end




    function out = isLinear(obj)
      out = true;
    end
  end

  methods  (Access = private)
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
      segmX = diff(obj.mesh.grid.X);
      segmY = diff(obj.mesh.grid.Y);
      segmZ = diff(obj.mesh.grid.Z);

      ncells = prod(obj.mesh.ncells);
      cellNeigh = zeros(ncells,6,'uint64');
      cellAreas = zeros(ncells,3);
      for k = 1:obj.mesh.ncells(3)
        for i = 1:obj.mesh.ncells(1)
          for j = 1:obj.mesh.ncells(2)
            cellId = obj.mesh.mapCellIds(i, j, k);
            % Neighbor Mapping:
            if i > 1
                cellNeigh(cellId, 3) = obj.mesh.mapCellIds(i-1,j,k);
            end
            if i < obj.mesh.ncells(1)
                cellNeigh(cellId, 4) = obj.mesh.mapCellIds(i+1,j,k);
            end
            if j > 1
                cellNeigh(cellId, 1) = obj.mesh.mapCellIds(i,j-1,k);
            end
            if j < obj.mesh.ncells(2)
                cellNeigh(cellId, 2) = obj.mesh.mapCellIds(i,j+1,k);
            end
            if k > 1
                cellNeigh(cellId, 5) = obj.mesh.mapCellIds(i,j,k-1);
            end
            if k < obj.mesh.ncells(3)
                cellNeigh(cellId, 6) = obj.mesh.mapCellIds(i,j,k+1);
            end
            
            % Face Area
            cellAreas(cellId, 1) = segmY(j) * segmZ(k);
            cellAreas(cellId, 2) = segmX(i) * segmZ(k);
            cellAreas(cellId, 3) = segmX(i) * segmY(j);
          end
        end
      end

      % Computing the permeability
      permCell = zeros(ncells,3);
      for mats=1:obj.nmat
        tmpMat=obj.materials.getMaterial(mats).PorousRock.getPermVoigt();
        permCell=permCell+tmpMat(1:3).*obj.mesh.matfrac(:,mats);
      end
      ActiveCells = obj.mesh.mapCellIds ~= 0;
      permCell=permCell(ActiveCells,:);

      % correcting matfrac, cleaning fractions for non-existent cells
      obj.mesh.matfrac(~ActiveCells(:),:)=0.;

      % temporay variables
      maxFaces = (obj.mesh.ncells(1)+1)*(obj.mesh.ncells(2)+1)*(obj.mesh.ncells(3)+1);
      tmpNeigh = zeros(maxFaces,2,'uint64');
      isActive = zeros(maxFaces,1,'logical');
      tmpTrans = zeros(maxFaces,1);
      tmpRhs = zeros(maxFaces,1);
      gamma = obj.materials.getFluid().getFluidSpecWeight();

      % Computing transmissibilities and RHS related to the gravity term
      count = 1;
      for cell=1:ncells
        [~,~,zNeiA] = obj.mesh.getIJKfromCellID(cell);
        zNeiA=obj.mesh.grid.centerZ(zNeiA);
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
              [~,~,zNeiB] = obj.mesh.getIJKfromCellID(neighbor);
              zNeiB=obj.mesh.grid.centerZ(zNeiB);
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

    function computeHPInitial(obj)
      ncells = obj.mesh.numberActiveCells;
      mu = obj.materials.getFluid().getDynViscosity();

      % Computing the initial H matrix
      tmpVec = (1/mu).*obj.trans;
      sumDiagTrans = accumarray(obj.facesNeigh(:),repmat(tmpVec,[2,1]),[ncells,1]);
      obj.H = sparse([obj.facesNeigh(:,1); obj.facesNeigh(:,2); (1:ncells)'],...
        [obj.facesNeigh(:,2); obj.facesNeigh(:,1); (1:ncells)'],...
        [-tmpVec; -tmpVec; sumDiagTrans], ncells, ncells);

      % Computing the initial P matrix
      maxcells = prod(obj.mesh.ncells);
      alpha=zeros(maxcells,1);
      poro=zeros(maxcells,1);
      for mats=1:obj.nmat
        tmpMat=obj.materials.getMaterial(mats).ConstLaw.getRockCompressibility();
        alpha=alpha+tmpMat.*obj.mesh.matfrac(:,mats);
        tmpMat=obj.materials.getMaterial(mats).PorousRock.getPorosity();
        poro=poro+tmpMat.*obj.mesh.matfrac(:,mats);
      end

      % Mapping and Reorder
      mapCellsOn = obj.mesh.mapCellIds ~= 0;
      map = obj.mesh.mapCellIds(mapCellsOn);
      alpha=alpha(mapCellsOn);
      alpha=alpha(map);
      poro=poro(mapCellsOn);
      poro=poro(map);
      volsCell = obj.mesh.computeVols();
      
      beta = obj.materials.getFluid().getFluidCompressibility();
      PVal = (alpha+beta*poro).*volsCell;
      obj.P = PVal.*speye(ncells);
    end

  end



  methods (Static)

    function out = getField()
      out = "pressure";
    end

    function out = isSymmetric()
      out = obj.Solver.isSymmetric();
    end

  end

end


