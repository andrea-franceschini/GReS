classdef domainManagerSed < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties
    % domain variables
    coordX (:,1)   % <- x array
    coordY (:,1)   % <- y array
    coordZ (:,1)   % <- z array
    bathy  (:,:)   % <- bathymetry of the bottom surface

    cm     (:,3)   % <- center of mass of the undeformed cell
    dl     (:,3)   % <- dimension of the undeformed cell    
    dzInc          % <- dz increment for grow z array

    % top surface information
    lastLay (:,:)  % <- last height
    top_dS
    top_e
    top_e0

    % dof manager control
    dof      (:,:,:)
    cellsByCol (:,:)
    ndofs

    initColumn (:,1)  % indication of the necessity of initialize the column

    % material variables
    matfrac (:,:)                % Material fractions per cell

    % interface variables
    interfX = struct('neigh',[],'area',[],'cm',[])
    interfY = struct('neigh',[],'area',[],'cm',[])
    interfZ = struct('neigh',[],'area',[],'cm',[])
  end

  methods (Access = public)
    function [obj, data] = domainManagerSed(data,nmat)
      default = struct('Grid',[],...
        'NewCellHeightControl',[],...
        'Bathymetry',missing);
      params = readInput(default,data);

      grid = domainManagerSed.getGrid(params.Grid);
      obj.coordX = grid.coordX;
      obj.coordY = grid.coordY;
      obj.coordZ = grid.coordZ;
      obj.dzInc = params.NewCellHeightControl;

      % Cell center coordinates
      [xc, yc, zc] = ndgrid((obj.coordX(1:end-1) + obj.coordX(2:end))/2,...
        (obj.coordY(1:end-1) + obj.coordY(2:end))/2,...
        (obj.coordZ(1:end-1) + obj.coordZ(2:end))/2);
      obj.cm = [xc(:) yc(:) zc(:)];

      % Cell volumes
      [dx, dy, dz] = ndgrid(diff(obj.coordX),diff(obj.coordY),diff(obj.coordZ));
      obj.dl = [dx(:) dy(:) dz(:)];

      map = [length(grid.coordX) length(grid.coordY) 2];

      % Defining the bathymetry for each column.
      if ~ismissing(params.Bathymetry)
        obj.bathy = domainManagerSed.getBathy(params.Bathymetry,map(1:2)-1);
      else
        obj.bathy = zeros(map(1:2)-1);
      end
      obj.lastLay = zeros(map(1:2)-1);

      % Defining the degree of freedom
      obj.ndofs = prod(map-1);
      obj.dof = reshape(1:obj.ndofs,map-1);
      obj.cellsByCol = ones(map(1:2)-1);
      obj.initColumn(1:obj.ndofs,1) = true;

      % % Defining the initial domain material mesh
      obj.matfrac = zeros(obj.ndofs,nmat);
      obj.matfrac(:,grid.mat) = 1.;

      % open space for sediments
      obj.top_dS = zeros(map(1:2)-1);
      obj.top_e = zeros(map(1:2)-1);
      obj.top_e0 = zeros(map(1:2)-1);
    end

    function growLayer(obj,map)
      nNewDofs = sum(map);
      newDofs = obj.ndofs+1:obj.ndofs+nNewDofs;
      ncells = size(obj.dof,1:3);

      % check if is necessary to add a new layer in z
      atTop = obj.cellsByCol == ncells(3);
      newlayer = any(and(atTop(:),map));
      if newlayer
        obj.dof(:,:,end+1) = zeros(ncells(1:2));
        tmpZ = obj.coordZ(end)+obj.dzInc;
        obj.coordZ(end+1) = tmpZ;
        ncells(3) = ncells(3)+1;
      end

      % update the dof for the newest cells.
      idI = repmat((1:ncells(1))',ncells(2),1);
      idJ = repelem((1:ncells(2))',ncells(1));
      pos = sub2ind(ncells,idI(map),idJ(map),obj.cellsByCol(map)+1);
      obj.dof(pos) = newDofs';
      obj.ndofs = obj.ndofs + nNewDofs;
      obj.cellsByCol(map) = obj.cellsByCol(map)+1;
    end

    function getConnectivity(obj,deform)
      ncells = size(obj.dof,1:3);

      % Make connections in the z direction
      tmp=sum(obj.cellsByCol-1,"all");
      obj.interfZ.neigh = zeros(tmp,2);
      obj.interfZ.area = zeros(tmp,1);
      obj.interfZ.cm = zeros(tmp,3);
      countZ = 1;
      for refJ=1:ncells(2)
        for refI=1:ncells(1)
          pillarCount = obj.cellsByCol(refI,refJ)-1;
          if pillarCount <= 1
            continue
          end
          interval = countZ:countZ+pillarCount-1;
          obj.interfZ.neigh(interval,:) = [reshape(...
            obj.dof(refI,refJ,1:pillarCount), [], 1),...
            reshape(obj.dof(refI,refJ,2:pillarCount+1), [], 1)];
          dx = obj.coordX(refI+1)-obj.coordX(refI);
          dy = obj.coordY(refJ+1)-obj.coordY(refJ);
          obj.interfZ.area(interval,1) = dx*dy;
          obj.interfZ.cm(interval,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
          obj.interfZ.cm(interval,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
          obj.interfZ.cm(interval,3) = (obj.coordZ(1:pillarCount)+obj.coordZ(2:pillarCount+1))/2.;
          countZ = countZ+pillarCount;
        end
      end

      % Make connections in the x and y direction
      obj.interfX.neigh = [];
      obj.interfX.area = [];
      obj.interfX.cm = [];

      obj.interfY.neigh = [];
      obj.interfY.area = [];
      obj.interfY.cm = [];

      % Allocating memory
      memMax = 0;
      for refJ=1:ncells(2)-1
        for refI=1:ncells(1)-1
          memMax=memMax+obj.cellsByCol(refI,refJ)+obj.cellsByCol(refI+1,refJ)-1;
          memMax=memMax+obj.cellsByCol(refI,refJ)+obj.cellsByCol(refI,refJ+1)-1;
        end
      end
      for refI=1:ncells(1)-1
        memMax=memMax+obj.cellsByCol(refI,ncells(2))+obj.cellsByCol(refI+1,ncells(2))-1;
      end
      for refJ=1:ncells(2)-1
        memMax=memMax+obj.cellsByCol(ncells(1),refJ)+obj.cellsByCol(ncells(1)-1,refJ+1)-1;
      end
      [zmin,zmax] = obj.getCoordinatesZ(deform);

      % Making connections
      neighX = zeros(memMax,2);
      areaX = zeros(memMax,1);
      cmfaceX = zeros(memMax,3);
      countX = 1;

      neighY = zeros(memMax,2);
      areaY = zeros(memMax,1);
      cmfaceY = zeros(memMax,3);
      countY = 1;

      if ncells(1)>1 && ncells(2)>1
        for refJ=1:ncells(2)-1
          for refI=1:ncells(1)-1
            dofA  = reshape(obj.dof(refI,refJ,1:obj.cellsByCol(refI,refJ)), [], 1);
            zMinA = zmin(dofA);
            zMaxA = zmax(dofA);

            % Finding conection in X direction.
            dofB  = reshape(obj.dof(refI+1,refJ,1:obj.cellsByCol(refI+1,refJ)), [], 1);
            zMinB = zmin(dofB);
            zMaxB = zmax(dofB);

            [neigh, face, faceCMz] = domainManagerSed.getPillarConnectivity(dofA, zMinA, zMaxA, dofB, zMinB, zMaxB);
            interval = reshape(countX:countX+length(face)-1, [], 1);
            neighX(interval,:) = neigh;
            areaX(interval,:) = face*(obj.coordY(refJ+1)-obj.coordY(refJ));
            cmfaceX(interval,1) = obj.coordX(refI+1);
            cmfaceX(interval,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
            cmfaceX(interval,3) = faceCMz;
            countX=countX+length(face);

            % Finding conection in Y direction.
            dofB  = reshape(obj.dof(refI,refJ+1,1:obj.cellsByCol(refI,refJ+1)), [], 1);
            zMinB = zmin(dofB);
            zMaxB = zmax(dofB);

            [neigh, face, faceCMz] = domainManagerSed.getPillarConnectivity(dofA, zMinA, zMaxA, dofB, zMinB, zMaxB);
            interval = reshape(countY:countY+length(face)-1, [], 1);
            neighY(interval,:) = neigh;
            areaY(interval,:) = face*(obj.coordX(refI+1)-obj.coordX(refI));
            cmfaceY(interval,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
            cmfaceY(interval,2) = obj.coordY(refJ+1);
            cmfaceY(interval,3) = faceCMz;
            countY=countY+length(face);
          end
        end
      end

      % Completing conection in X direction.
      if ncells(1)>1
        for refI=1:ncells(1)-1
          dofA  = reshape(obj.dof(refI,end,1:obj.cellsByCol(refI,end)), [], 1);
          zMinA = zmin(dofA);
          zMaxA = zmax(dofA);

          dofB  = reshape(obj.dof(refI+1,end,1:obj.cellsByCol(refI+1,end)), [], 1);
          zMinB = zmin(dofB);
          zMaxB = zmax(dofB);

          [neigh, face, faceCMz] = domainManagerSed.getPillarConnectivity(dofA, zMinA, zMaxA, dofB, zMinB, zMaxB);
          interval = reshape(countX:countX+length(face)-1, [], 1);
          neighX(interval,:) = neigh;
          areaX(interval,:) = face*(obj.coordY(end)-obj.coordY(end-1));
          cmfaceX(interval,1) = obj.coordX(refI+1);
          cmfaceX(interval,2) = (obj.coordY(end)+obj.coordY(end-1))/2.;
          cmfaceX(interval,3) = faceCMz;
          countX=countX+length(face);
        end
      end

      % Remove unused preallocated entries
      obj.interfX.neigh = neighX(1:countX-1,:);
      obj.interfX.area = areaX(1:countX-1,:);
      obj.interfX.cm = cmfaceX(1:countX-1,:);

      % Finding conection in Y direction.
      if ncells(2)>1
        for refJ=1:ncells(2)-1
          dofA  = reshape(obj.dof(end,refJ,1:obj.cellsByCol(end,refJ)), [], 1);
          zMinA = zmin(dofA);
          zMaxA = zmax(dofA);

          dofB  = reshape(obj.dof(end,refJ+1,1:obj.cellsByCol(end,refJ+1)), [], 1);
          zMinB = zmin(dofB);
          zMaxB = zmax(dofB);

          [neigh, face, faceCMz] = domainManagerSed.getPillarConnectivity(dofA, zMinA, zMaxA, dofB, zMinB, zMaxB);
          interval = reshape(countY:countY+length(face)-1, [], 1);
          neighY(interval,:) = neigh;
          areaY(interval,:) = face*(obj.coordX(end)-obj.coordX(end-1));
          cmfaceY(interval,1) = (obj.coordX(end)+obj.coordX(end-1))/2.;
          cmfaceY(interval,2) = obj.coordY(refJ+1);
          cmfaceY(interval,3) = faceCMz;
          countY=countY+length(face);
        end
      end

      % Remove unused preallocated entries
      obj.interfY.neigh = neighY(1:countY-1,:);
      obj.interfY.area = areaY(1:countY-1,:);
      obj.interfY.cm = cmfaceY(1:countY-1,:);
    end

    function getConnectivityPARALLEL(obj,deform)
      ncells = size(obj.dof,1:3);
      nx = ncells(1);
      ny = ncells(2);

      dof_pt = obj.dof;
      dofByCol_pt = obj.cellsByCol;

      coordX_pt = obj.coordX;
      coordY_pt = obj.coordY;
      coordZ_pt = obj.coordZ;

      [zmin,zmax] = obj.getCoordinatesZ(deform);

      % Z-direction connectivity
      nPillars = nx*ny;
      neighCell = cell(nPillars,1);
      areaCell  = cell(nPillars,1);
      cmCell    = cell(nPillars,1);

      parfor p = 1:nPillars
        [refI,refJ] = ind2sub([nx,ny],p);
        pillarCount = dofByCol_pt(refI,refJ)-1;
        if pillarCount < 1
          neighCell{p} = zeros(0,2);
          areaCell{p}  = zeros(0,1);
          cmCell{p}    = zeros(0,3);
          continue
        end

        dofs = reshape( ...
          dof_pt(refI,refJ,1:pillarCount+1),[],1);
        neighCell{p} = [dofs(1:end-1),dofs(2:end)];

        dx = coordX_pt(refI+1)-coordX_pt(refI);
        dy = coordY_pt(refJ+1)-coordY_pt(refJ);
        areaCell{p} = repmat(dx*dy,pillarCount,1);
        cmCell{p} = [ ...
          repmat((coordX_pt(refI+1)+coordX_pt(refI))/2,pillarCount,1), ...
          repmat((coordY_pt(refJ+1)+coordY_pt(refJ))/2,pillarCount,1), ...
          (coordZ_pt(1:pillarCount)+coordZ_pt(2:pillarCount+1))/2 ...
          ];
      end

      obj.interfZ.neigh = vertcat(neighCell{:});
      obj.interfZ.area  = vertcat(areaCell{:});
      obj.interfZ.cm    = vertcat(cmCell{:});

      % X-direction connectivity

      nPairsX = (nx-1)*ny;

      neighCell = cell(nPairsX,1);
      areaCell  = cell(nPairsX,1);
      cmCell    = cell(nPairsX,1);

      parfor p = 1:nPairsX

        [refI,refJ] = ind2sub([nx-1,ny],p);

        nA = dofByCol_pt(refI,refJ);
        nB = dofByCol_pt(refI+1,refJ);

        dofA = reshape(dof_pt(refI,refJ,1:nA),[],1);
        dofB = reshape(dof_pt(refI+1,refJ,1:nB),[],1);

        [neigh,face,faceCMz] = ...
          domainManagerSed.getPillarConnectivity( ...
          dofA,zmin(dofA),zmax(dofA), ...
          dofB,zmin(dofB),zmax(dofB));

        nFaces = numel(face);

        dy = coordY_pt(refJ+1)-coordY_pt(refJ);

        neighCell{p} = neigh;
        areaCell{p}  = face*dy;

        cmCell{p} = [ ...
          repmat(coordX_pt(refI+1),nFaces,1), ...
          repmat((coordY_pt(refJ+1)+coordY_pt(refJ))/2,nFaces,1), ...
          faceCMz ...
          ];

      end

      % Preserve the original ordering: refJ outer, refI inner
      neighCell = reshape(neighCell,nx-1,ny).';
      areaCell  = reshape(areaCell,nx-1,ny).';
      cmCell    = reshape(cmCell,nx-1,ny).';

      obj.interfX.neigh = vertcat(neighCell{:});
      obj.interfX.area  = vertcat(areaCell{:});
      obj.interfX.cm    = vertcat(cmCell{:});

      % Y-direction connectivity

      nPairsY = nx*(ny-1);

      neighCell = cell(nPairsY,1);
      areaCell  = cell(nPairsY,1);
      cmCell    = cell(nPairsY,1);

      parfor p = 1:nPairsY

        [refI,refJ] = ind2sub([nx,ny-1],p);

        nA = dofByCol_pt(refI,refJ);
        nB = dofByCol_pt(refI,refJ+1);

        dofA = reshape(dof_pt(refI,refJ,1:nA),[],1);
        dofB = reshape(dof_pt(refI,refJ+1,1:nB),[],1);

        [neigh,face,faceCMz] = ...
          domainManagerSed.getPillarConnectivity( ...
          dofA,zmin(dofA),zmax(dofA), ...
          dofB,zmin(dofB),zmax(dofB));

        nFaces = numel(face);

        dx = coordX_pt(refI+1)-coordX_pt(refI);

        neighCell{p} = neigh;
        areaCell{p}  = face*dx;

        cmCell{p} = [ ...
          repmat((coordX_pt(refI+1)+coordX_pt(refI))/2,nFaces,1), ...
          repmat(coordY_pt(refJ+1),nFaces,1), ...
          faceCMz ...
          ];

      end

      % Preserve the original ordering: refJ outer, refI inner
      neighCell = reshape(neighCell,nx,ny-1).';
      areaCell  = reshape(areaCell,nx,ny-1).';
      cmCell    = reshape(cmCell,nx,ny-1).';

      obj.interfY.neigh = vertcat(neighCell{:});
      obj.interfY.area  = vertcat(areaCell{:});
      obj.interfY.cm    = vertcat(cmCell{:});
    end
    
 
    function out = getMapFromDofs(obj,dofs)
      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      map = ismember(obj.dof,dofs);
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [idI,idJ,~]=ind2sub(size(obj.dof),pos(id));
      out = sub2ind(size(obj.bathy), idI', idJ');
    end

    function [dofs] = getBord(obj,label)
      % GETBORDDOFS Returns dofs for the border.
      %
      % Outputs:
      %   id        - DOFs at the boundary
      %
      % Supported labels:
      %   'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'

      switch lower(label)
        case "xmin"
          tmp = obj.dof(1,:,:);
          dofs(:) = tmp(tmp~=0);
        case "xmax"
          tmp = obj.dof(end,:,:);
          dofs(:) = tmp(tmp~=0);
        case "ymin"
          tmp = obj.dof(:,1,:);
          dofs(:) = tmp(tmp~=0);
        case "ymax"
          tmp = obj.dof(:,end,:);
          dofs(:) = tmp(tmp~=0);
        case "zmin"
          tmp = obj.dof(:,:,1);
          dofs(:) = tmp(tmp~=0);
        case "zmax"
          ncells = size(obj.dof,1:3);
          idI = repmat((1:ncells(1))',ncells(2),1);
          idJ = repelem((1:ncells(2))',ncells(1));
          pos = sub2ind(ncells,idI,idJ,obj.cellsByCol(:));
          tmp = obj.dof(pos);
          dofs(:) = tmp(tmp~=0);
        otherwise
          dofs = [];
      end
      dofs=dofs';
    end

    function updateConectivities(obj,map,dz)
      ncells = size(obj.dof,1:3);
      nNewDofs = sum(map);
      obj.growLayer(map);
      obj.lastLay(map) = dz;

      % Make connections in the z direction
      map=reshape(map,ncells(1:2));
      neigh = zeros(nNewDofs,2);
      area = zeros(nNewDofs,1);
      cm = zeros(nNewDofs,3);
      countZ = 1;
      for refJ=1:ncells(2)
        for refI=1:ncells(1)
          if map(refI,refJ)
            refK = obj.cellsByCol(refI,refJ);
            neigh(countZ,:) = [obj.dof(refI,refJ,refK) obj.dof(refI,refJ,refK)];

            dx = obj.coordX(refI+1)-obj.coordX(refI);
            dy = obj.coordY(refJ+1)-obj.coordY(refJ);
            area(countZ,1) = dx*dy;
            cm(countZ,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
            cm(countZ,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
            cm(countZ,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
            countZ = countZ+1;
          end
        end
      end
      obj.interfZ.neigh = [obj.interfZ.neigh; neigh];
      obj.interfZ.area = [obj.interfZ.area; area];
      obj.interfZ.cm = [obj.interfZ.cm; cm];

      % Make connections in the x and y direction
      countX = 1;
      neighX = zeros((ncells(1)-1)*ncells(2)*2,2);
      areaX = zeros(nNewDofs,1);
      cmX = zeros(nNewDofs,3);
      countY = 1;
      neighY = zeros((ncells(2)-1)*ncells(1)*2,2);
      areaY = zeros(nNewDofs,1);
      cmY = zeros(nNewDofs,3);
      for refJ=1:ncells(2)
        dy = obj.coordY(refJ+1)-obj.coordY(refJ);
        for refI=1:ncells(1)
          dx = obj.coordX(refI+1)-obj.coordX(refI);
          dz = obj.coordZ(refK+1)-obj.coordZ(refK);
          if map(refI,refJ)
            refK = obj.cellsByCol(refI,refJ);
            dofRef = obj.dof(refI,refJ,refK);

            if refI==1
              dofA = obj.dof(refI+1,refJ,refK);
              if dofA ~= 0
                neighX(countX,:) = [dofRef dofA];
                areaX(countX,1) = dy*dz;
                cmX(countX,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmX(countX,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmX(countX,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countX=countX+1;
              end
            elseif refI==ncells(1)
              dofA = obj.dof(refI-1,refJ,refK);
              if ~ismember([dofA dofRef], obj.interfX.neigh, 'rows') && ...
                  ~ismember([dofA dofRef], neighX, 'rows') && dofA ~=0
                neighX(countX,:) = [dofA dofRef];
                areaX(countX,1) = dy*dz;
                cmX(countX,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmX(countX,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmX(countX,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countX=countX+1;
              end
            else
              dofA = obj.dof(refI-1,refJ,refK);
              if ~ismember([dofA dofRef], obj.interfX.neigh, 'rows') && ...
                  ~ismember([dofA dofRef], neighX, 'rows') && dofA ~=0
                neighX(countX,:) = [dofA dofRef];
                areaX(countX,1) = dy*dz;
                cmX(countX,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmX(countX,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmX(countX,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countX=countX+1;
              end
              dofA = obj.dof(refI+1,refJ,refK);
              if dofA ~=0
                neighX(countX,:) = [dofRef dofA];
                areaX(countX,1) = dy*dz;
                cmX(countX,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmX(countX,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmX(countX,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countX=countX+1;
              end
            end

            if refJ==1
              dofA = obj.dof(refI,refJ+1,refK);
              if dofA ~=0
                neighY(countY,:) = [dofRef dofA];
                areaY(countY,1) = dx*dz;
                cmY(countY,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmY(countY,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmY(countY,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countY=countY+1;
              end
            elseif refJ==ncells(2)
              dofA = obj.dof(refI,refJ-1,refK);
              if ~ismember([dofA dofRef], obj.interfY.neigh, 'rows') && ...
                  ~ismember([dofA dofRef], neighY, 'rows') && dofA ~=0
                neighY(countY,:) = [dofA dofRef];
                areaY(countY,1) = dx*dz;
                cmY(countY,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmY(countY,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmY(countY,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countY=countY+1;
              end
            else
              dofA = obj.dof(refI,refJ-1,refK);
              if ~ismember([dofA dofRef], obj.interfY.neigh, 'rows') && ...
                  ~ismember([dofA dofRef], neighY, 'rows') && dofA ~=0
                neighY(countY,:) = [dofA dofRef];
                areaY(countY,1) = dx*dz;
                cmY(countY,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmY(countY,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmY(countY,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countY=countY+1;
              end
              dofA = obj.dof(refI,refJ+1,refK);
              if dofA ~=0
                neighY(countY,:) = [dofRef dofA];
                areaY(countY,1) = dx*dz;
                cmY(countY,1) = (obj.coordX(refI+1)+obj.coordX(refI))/2.;
                cmY(countY,2) = (obj.coordY(refJ+1)+obj.coordY(refJ))/2.;
                cmY(countY,3) = (obj.coordZ(refK)+obj.coordZ(refK))/2.;
                countY=countY+1;
              end
            end
          end
        end
      end

      obj.interfX.neigh(end+1:end+countX-1,:)=neighX(1:countX-1,:);
      obj.interfX.area(end+1:end+countX-1,:)=areaX(1:countX-1,:);
      obj.interfX.cm(end+1:end+countX-1,:)=cmX(1:countX-1,:);
      obj.interfY.neigh(end+1:end+countY-1,:)=neighY(1:countY-1,:);
      obj.interfY.area(end+1:end+countY-1,:)=areaY(1:countY-1,:);
      obj.interfY.cm(end+1:end+countY-1,:)=cmY(1:countY-1,:);
    end

    function vol = getCellVols(obj,dl,dofs)
      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      vol = prod(obj.dl(dofs,1:2),2).*(obj.dl(dofs,3)+dl(dofs));
    end

    function [trans,neigh] = getFaceTransmissibility(obj,perm,deform)

      getConnectivity(obj,deform);
      defCM = getDeformCM(obj,deform);
      ncells(1) = length(obj.interfX.area);
      ncells(2) = length(obj.interfY.area);
      ncells(3) = length(obj.interfZ.area);
      trans(1:sum(ncells),1)=0;
      neigh(1:sum(ncells),2)=0;

      if ~isempty(obj.interfX.area)
        interval = 1:ncells(1);
        trans(interval) = domainManagerSed.halfTransmissibility(obj.interfX,...
          [1 0 0],perm,defCM);
        neigh(interval,:) = obj.interfX.neigh;
      end

      if ~isempty(obj.interfY.area)
        interval = ncells(1)+1:sum(ncells(1:2));
        trans(interval) = domainManagerSed.halfTransmissibility(obj.interfY,...
          [0 1 0],perm,defCM);
        neigh(interval,:) = obj.interfY.neigh;
      end

      if ~isempty(obj.interfZ.area)
        interval = sum(ncells(1:2))+1:sum(ncells);
        trans(interval) = domainManagerSed.halfTransmissibility(obj.interfZ,...
          [0 0 1],perm,defCM);
        neigh(interval,:) = obj.interfZ.neigh;
      end

    end

    function [trans,neigh] = getFaceTransmissibilityPARALLEL(obj,perm,deform)
      getConnectivity(obj,deform);
      defCM = getDeformCM(obj,deform);

      interfaces = {obj.interfX, obj.interfY, obj.interfZ};
      normals = [1 0 0; 0 1 0; 0 0 1];

      transCell = cell(3,1);
      neighCell = cell(3,1);

      parfor i = 1:3
        interface = interfaces{i};
        if ~isempty(interface.area)
          transCell{i} = domainManagerSed.halfTransmissibility( ...
            interface, normals(i,:), perm, defCM);
          neighCell{i} = interface.neigh;
        else
          transCell{i} = zeros(0,1);
          neighCell{i} = zeros(0,2);
        end
      end

      trans = vertcat(transCell{:});
      neigh = vertcat(neighCell{:});
    end

    function trans = getBcHalfTrans(obj,type,perm,deform,dofs)
      % COMPUTEHALFTRANS Computes cells half-transmissibility.

      % Get mesh dimension and update with the deformation
      dx = obj.dl(dofs,1);
      dy = obj.dl(dofs,2);
      dz = obj.dl(dofs,3) + deform;

      % Computing half transmissibilities
      switch type
        case 1
          trans = (dy.*dz)./(dx/2).*perm(:,1);      
        case 2
          trans = (dx.*dz)./(dy/2).*perm(:,2);      
        case 3
          trans = (dx.*dy)./(dz/2).*perm(:,3);
      end
    end
  end

  methods (Access = private)
    function [zmin,zmax] = getCoordinatesZ(obj,deform)
      ncells = size(obj.dof,1:3);

      zmin(1:obj.ndofs,1) = 0.;
      zmax(1:obj.ndofs,1) = 0.;
      for refJ=1:ncells(2)
        for refI=1:ncells(1)
          column = reshape(obj.dof(refI,refJ,:), [], 1);
          dofs = column(column~=0);
          ndofsCol = nnz(column);
          if ndofsCol==1
            posRef = obj.bathy(refI,refJ)+deform(dofs);
            zmin(dofs,1) = posRef+obj.coordZ(1:ndofsCol,1);
            zmax(dofs,1) = posRef+obj.coordZ(2:ndofsCol+1,1);
          else
            posRef = obj.bathy(refI,refJ)+deform(dofs);
            zmin(dofs(1:end-1),1) = posRef(1:end-1)+obj.coordZ(1:ndofsCol-1,1);
            zmax(dofs(1:end-1),1) = posRef(1:end-1)+obj.coordZ(2:ndofsCol,1);

            zmin(dofs(end),1) = posRef(end)+obj.coordZ(ndofsCol,1);
            zmax(dofs(end),1) = posRef(end)+obj.coordZ(ndofsCol,1)+obj.lastLay(refI,refJ);
          end
        end
      end
    end

    function defCM = getDeformCM(obj,deform)
      % The deformation of the last cell is always zero, because is a
      % special cell that grow in function of the material being deposited.
      defCM = obj.cm;
      ncells = size(obj.bathy);
      for j=1:ncells(2)
        for i=1:ncells(1)
          dofs = obj.dof(i,j,:);
          dofs = dofs(dofs~=0);
          dofs = dofs(1:end);
          if ~isempty(dofs)
            acc = 0;
            for k=1:length(dofs)
              dofRef=dofs(k);
              defCM(dofRef,3)=defCM(dofRef,3)+deform(dofRef)/2.+acc;
              acc=acc+deform(dofRef);
            end            
          end
        end
      end
    end

  end

  methods (Static)
    function data = getGrid(input)
      % Contruct grid.
      default = struct('type',string.empty);
      params = readInput(default,input);

      switch lower(params.type)
        case "classic"
          default = struct('division',double.empty,...
            'size',double.empty,'materialTag',double.empty);
          params = readInput(default,params);

          ncells = params.division;
          dim = params.size;
          data.coordX = linspace(0,dim(1),ncells(1)+1);
          data.coordY = linspace(0,dim(2),ncells(2)+1);
          data.coordZ = linspace(0,dim(3),2);

          data.mat = params.materialTag;
        case "explicit"
          default = struct('xcoord',missing,...
            'xfile',missing,...
            'ycoord',missing,...
            'yfile',missing,...
            'zsize',double.empty,...
            'materialTag',double.empty);
          params = readInput(default,params.Grid);
          data.mat = params.materialTag;

          % x coordanate
          if ~ismissing(params.xcoord) && ~ismissing(params.xfile)
            gresLog().error("The grid parameters for the simulation is not well defined!");
          elseif ~ismissing(xfile)
            data.coordX = load(params.xfile);
          else
            data.coordX = params.xcoord;
          end

          % y coordanate
          if ~ismissing(params.ycoord) && ~ismissing(params.yfile)
            gresLog().error("The grid parameters for the simulation is not well defined!");
          elseif ~ismissing(yfile)
            data.coordY = load(params.yfile);
          else
            data.coordY = params.ycoord;
          end

          % z coordanate          
          if ~ismissing(params.zsize)
            gresLog().error("The grid parameters for the simulation is not well defined!");
          else
            data.coordZ = linspace(0,params.zsize,2);
          end
        otherwise
          gresLog().error("The grid parameters for the simulation is not well defined!");
      end
    end

    function data = getBathy(input,map)
      if ~isfile(input.file)
        gresLog().error("Bathymetry file not found!");
      end
      vls = readmatrix(input.file);

      % Target dimensions
      nx = map(1);
      ny = map(2);

      % Input dimensions
      [mx,my] = size(vls);

      % Same dimensions
      if mx == nx && my == ny
        data = vls;
        return;
      end

      % Single value
      if isscalar(vls)
        data = repmat(vls,nx,ny);
        return;
      end

      % Normalized coordinates
      xnew = linspace(0,1,nx);
      ynew = linspace(0,1,ny);

      % One-dimensional input
      if isvector(vls)
        n = numel(vls);
        xref = linspace(0,1,n);

        if n == 2
          polyType = 'linear';
        elseif n == 3
          polyType = 'cubic';
        else
          polyType = 'spline';
        end

        if mx > 1
          % Interpolation along first dimension
          tmp = interp1(xref,vls(:),xnew,polyType);
          data = repmat(tmp(:),1,ny);
        else
          % Interpolation along second dimension
          tmp = interp1(xref,vls(:),ynew,polyType);
          data = repmat(tmp(:)',nx,1);
        end
        return;
      end

      % Two-dimensional input
      if min(mx,my) == 2
        polyType = 'linear';
      elseif min(mx,my) == 3
        polyType = 'cubic';
      else
        polyType = 'spline';
      end
      xref = linspace(0,1,mx);
      yref = linspace(0,1,my);
      [Xref,Yref] = ndgrid(xref,yref);
      [Xnew,Ynew] = ndgrid(xnew,ynew);
      data = interp2(Yref,Xref,vls,Ynew,Xnew,polyType);
    end

    function [neigh, faceDz, faceCMz] = getPillarConnectivity(dofA, zMinA, zMaxA,...
        dofB, zMinB, zMaxB)
      %GETPILLARCONNECTIVITY Compute lateral connectivity between two pillars.
      %
      %   [neigh, faceLength, faceCenterZ] = getPillarConnectivity( ...
      %       dofA, zMinA, zMaxA, dofB, zMinB, zMaxB)
      %
      %   Uses a two-pointer algorithm to identify overlapping vertical
      %   intervals between cells belonging to two adjacent pillars.
      %
      %   INPUTS:
      %     dofA, dofB   - Global cell DOF indices for each pillar.
      %     zMinA, zMinB - Lower elevations of the cells.
      %     zMaxA, zMaxB - Upper elevations of the cells.
      %
      %   OUTPUTS:
      %     neigh       - [N x 2] global DOF connectivity.
      %     faceLength  - [N x 1] vertical overlap lengths.
      %     faceCenterZ - [N x 1] elevations of contact-face centers.
      %
      %   NOTES:
      %     Cells must be sorted from bottom to top within each pillar.
      %     Cell intervals must not overlap within the same pillar.
      %     Connections with zero overlap are excluded.

      % Ensure column vectors
      dofA = dofA(:);
      dofB = dofB(:);
      zMinA = zMinA(:);
      zMaxA = zMaxA(:);
      zMinB = zMinB(:);
      zMaxB = zMaxB(:);

      nA = numel(dofA);
      nB = numel(dofB);

      % Handle empty pillars
      if nA == 0 || nB == 0
        neigh = zeros(0, 2);
        faceDz = zeros(0, 1);
        faceCMz = zeros(0, 1);
        return
      end

      % Preallocate maximum possible connections
      nMax = nA + nB - 1;

      neigh = zeros(nMax, 2);
      faceDz = zeros(nMax, 1);
      faceCMz = zeros(nMax, 1);

      % Initialize pointers and connection counter
      a = 1;
      b = 1;
      n = 0;

      while a <= nA && b <= nB

        % Compute overlapping interval
        lower = max(zMinA(a), zMinB(b));
        upper = min(zMaxA(a), zMaxB(b));

        overlap = upper - lower;

        if overlap > 0
          n = n + 1;

          neigh(n,:) = [dofA(a), dofB(b)];
          faceDz(n) = overlap;
          faceCMz(n) = (lower + upper)/2;
        end

        % Advance the pointer of the cell ending first
        if zMaxA(a) < zMaxB(b)
          a = a + 1;
        elseif zMaxA(a) > zMaxB(b)
          b = b + 1;
        else
          a = a + 1;
          b = b + 1;
        end
      end

      % Remove unused preallocated entries
      neigh = neigh(1:n,:);
      faceDz = faceDz(1:n);
      faceCMz = faceCMz(1:n);
    end

    function trans = halfTransmissibility(interface,normal,perm,defCM)
      Cik = interface.cm-defCM(interface.neigh(:,1),:);
      Nik = normal.*interface.area;
      Ki = perm(interface.neigh(:,1),:);
      Tik = sum(Nik .* Ki .* Cik, 2) ./ sum(Cik.^2, 2);

      Cik = interface.cm-defCM(interface.neigh(:,2),:);
      Nik = -normal.*interface.area;
      Ki = perm(interface.neigh(:,2),:);
      Tki = sum(Nik .* Ki .* Cik, 2) ./ sum(Cik.^2, 2);
      trans = 1./(1./Tik+1./Tki);
    end



  end

end