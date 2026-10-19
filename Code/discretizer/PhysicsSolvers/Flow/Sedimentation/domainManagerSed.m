classdef domainManagerSed < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties
    % domain variables
    coordX (:,1)   % <- x array
    coordY (:,1)   % <- y array
    coordZ (:,1)   % <- z array
    bathy  (:,:)   % <- bathymetry of the bottom surface
    dzInc          % <- dz increment for grow z array
    lastLay (:,:)  % <- last height

    % dof manager control
    dof      (:,:,:)
    dofByCol (:,:)
    ndofs

    % interface variables
    interfX = struct('neigh',[],'area',[],'cm',[])
    interfY = struct('neigh',[],'area',[],'cm',[])
    interfZ = struct('neigh',[],'area',[],'cm',[])
  end

  methods (Access = public)
    function obj = domainManagerSed(data)
      default = struct('Grid',[],...
        'Initial',[],...
        'NewCellHeightControl',[],...
        'bathy',missing);
      params = readInput(default,data);

      grid = domainManagerSed.getGrid(params.Grid);
      obj.coordX = grid.coordX;
      obj.coordY = grid.coordY;
      obj.coordZ = grid.coordZ;
      obj.dzInc = params.NewCellHeightControl;

      map = [length(grid.coordX) length(grid.coordY) length(grid.coordZ)];

      % Defining the bathymetry for each column.
      if ~ismissing(params.bathy)
        obj.bathy = domainManagerSed.getBathy(params.bathy,map(1:2)-1);
      else
        obj.bathy = zeros([length(obj.coordX) length(obj.coordY)]);
      end
      obj.lastLay = zeros([length(obj.coordX)-1 length(obj.coordY)-1]);

      % Defining the degree of freedom
      obj.dof = domainManagerSed.getInitialDofs(params.Initial,map-1);
      obj.ndofs = nnz(obj.dof);
      obj.dofByCol = zeros(map(1:2)-1);
      if map(3)==1
        obj.dofByCol(:) = 1;
      else
        for refI=1:map(1)-1
          for refJ=1:map(2)-1
            obj.dofByCol(refI,refJ) = nnz(obj.dof(refI,refJ,:));
          end
        end
      end
    end

    function growLayer(obj,map)
      nNewDofs = sum(map);
      newDofs = obj.ndofs+1:obj.ndofs+nNewDofs;
      ncells = size(obj.dof,1:3);

      % check if is necessary to add a new layer in z
      atTop = obj.dofByCol == ncells(3);
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
      pos = sub2ind(ncells,idI(map),idJ(map),obj.dofByCol(map)+1);
      obj.dof(pos) = newDofs';
      obj.ndofs = obj.ndofs + nNewDofs;
      obj.dofByCol(map) = obj.dofByCol(map)+1;
    end

    function getConnectivity(obj,deform)
      ncells = size(obj.dof,1:3);

      % Make connections in the z direction
      tmp=sum(obj.dofByCol-1,"all");
      obj.interfZ.neigh = zeros(tmp,2);
      obj.interfZ.area = zeros(tmp,1);
      obj.interfZ.cm = zeros(tmp,3);
      countZ = 1;
      for refJ=1:ncells(2)
        for refI=1:ncells(1)
          pillarCount = obj.dofByCol(refI,refJ)-1;
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
          memMax=memMax+obj.dofByCol(refI,refJ)+obj.dofByCol(refI+1,refJ)-1;
          memMax=memMax+obj.dofByCol(refI,refJ)+obj.dofByCol(refI,refJ+1)-1;
        end
      end
      for refI=1:ncells(1)-1
        memMax=memMax+obj.dofByCol(refI,ncells(2))+obj.dofByCol(refI+1,ncells(2))-1;
      end
      for refJ=1:ncells(2)-1
        memMax=memMax+obj.dofByCol(ncells(1),refJ)+obj.dofByCol(ncells(1)-1,refJ+1)-1;
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
            dofA  = reshape(obj.dof(refI,refJ,1:obj.dofByCol(refI,refJ)), [], 1);
            zMinA = zmin(dofA);
            zMaxA = zmax(dofA);

            % Finding conection in X direction.
            dofB  = reshape(obj.dof(refI+1,refJ,1:obj.dofByCol(refI+1,refJ)), [], 1);
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
            dofB  = reshape(obj.dof(refI,refJ+1,1:obj.dofByCol(refI,refJ+1)), [], 1);
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
          dofA  = reshape(obj.dof(refI,end,1:obj.dofByCol(refI,end)), [], 1);
          zMinA = zmin(dofA);
          zMaxA = zmax(dofA);

          dofB  = reshape(obj.dof(refI+1,end,1:obj.dofByCol(refI+1,end)), [], 1);
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
          dofA  = reshape(obj.dof(end,refJ,1:obj.dofByCol(end,refJ)), [], 1);
          zMinA = zmin(dofA);
          zMaxA = zmax(dofA);

          dofB  = reshape(obj.dof(end,refJ+1,1:obj.dofByCol(end,refJ+1)), [], 1);
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

    % % % % function makeMesh(obj,a)
    % % % % end

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
          pos = sub2ind(ncells,idI,idJ,obj.dofByCol(:));
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
            refK = obj.dofByCol(refI,refJ);
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
            refK = obj.dofByCol(refI,refJ);
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

    % % % % function [coord,conect] = getMesh1(obj,a)
    % % % %   npoint = sum(4*idxTop+4,"all");
    % % % %   ncell = sum(idxTop,"all");
    % % % %
    % % % %   coord(1:npoint,1:3) = 0.;
    % % % %   conec(1:ncell,1:8)  = 0;
    % % % %
    % % % %   count = 0;
    % % % %   for col=1:prod(obj.ncells(1:2))
    % % % %     idI = mod(col-1,obj.ncells(1))+1;
    % % % %     idJ = (col-idI)/obj.ncells(1)+1;
    % % % %
    % % % %     coord(4*count+1,1) = coordX(idI);
    % % % %     coord(4*count+2,1) = coordX(idI+1);
    % % % %     coord(4*count+3,1) = coordX(idI+1);
    % % % %     coord(4*count+4,1) = coordX(idI);
    % % % %
    % % % %     coord(4*count+1,2) = coordY(idJ);
    % % % %     coord(4*count+2,2) = coordY(idJ);
    % % % %     coord(4*count+3,2) = coordY(idJ+1);
    % % % %     coord(4*count+4,2) = coordY(idJ+1);
    % % % %
    % % % %     coord(4*count+1:4*count+4,3) = coordZ(1);
    % % % %     count = count+1;
    % % % %     topLay = idxTop(col);
    % % % %     for lay = 1:topLay
    % % % %       coord(4*count+1,1) = coordX(idI);
    % % % %       coord(4*count+2,1) = coordX(idI+1);
    % % % %       coord(4*count+3,1) = coordX(idI+1);
    % % % %       coord(4*count+4,1) = coordX(idI);
    % % % %
    % % % %       coord(4*count+1,2) = coordY(idJ);
    % % % %       coord(4*count+2,2) = coordY(idJ);
    % % % %       coord(4*count+3,2) = coordY(idJ+1);
    % % % %       coord(4*count+4,2) = coordY(idJ+1);
    % % % %
    % % % %       if lay == topLay
    % % % %         if ~initCol(col)
    % % % %           coord(4*count+1:4*count+4,3) = coordZ(lay)+sedAcc(col);
    % % % %         else
    % % % %           coord(4*count+1:4*count+4,3) = coordZ(lay+1);
    % % % %         end
    % % % %       else
    % % % %         coord(4*count+1:4*count+4,3) = coordZ(lay+1);
    % % % %       end
    % % % %       cell = dofs(obj.dof(idI,idJ,lay));
    % % % %       conec(cell,1)=4*(count-1)+1;
    % % % %       conec(cell,2)=4*(count-1)+2;
    % % % %       conec(cell,3)=4*(count-1)+3;
    % % % %       conec(cell,4)=4*(count-1)+4;
    % % % %       conec(cell,5)=4*count+1;
    % % % %       conec(cell,6)=4*count+2;
    % % % %       conec(cell,7)=4*count+3;
    % % % %       conec(cell,8)=4*count+4;
    % % % %       count=count+1;
    % % % %     end
    % % % %   end
    % % % % end




  end

  methods (Static)
    function data = getGrid(input)
      % Contruct grid.
      default = struct('type',string.empty);
      params = readInput(default,input);

      switch lower(params.type)
        case "classic"
          default = struct('division',double.empty,...
            'size',double.empty);
          params = readInput(default,params);

          ncells = params.division;
          dim = params.size;
          data.coordX = linspace(0,dim(1),ncells(1)+1);
          data.coordY = linspace(0,dim(2),ncells(2)+1);
          data.coordZ = linspace(0,dim(3),ncells(3)+1);
        case "explicit"
          default = struct('xcoord',missing,...
            'xfile',missing,...
            'ycoord',missing,...
            'yfile',missing,...
            'zcoord',missing,...
            'zfile',missing);
          params = readInput(default,params.Grid);

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
          if ~ismissing(params.zcoord) && ~ismissing(params.zfile)
            gresLog().error("The grid parameters for the simulation is not well defined!");
          elseif ~ismissing(zfile)
            data.coordZ = load(params.zfile);
          else
            data.coordZ = params.zcoord;
          end
        otherwise
          gresLog().error("The grid parameters for the simulation is not well defined!");
      end
    end

    function data = getBathy(input,map)
      vls = readmatrix(input.file,'Delimiter','\t');
      mapDim = size(vls);
      if all(mapDim==map)
        tmp = vls;
      else
        [x_ref, y_ref] = meshgrid(0:1/(mapDim(1)-1):1, 0:1/(mapDim(2)-1):1);
        [x_new, y_new] = meshgrid(0:1/(map(1)-1):1, 0:1/(map(2)-1):1);
        tmp = interp2(x_ref, y_ref, vls, x_new, y_new, 'cubic')';
      end
      data = tmp;
    end

    function data = getInitialDofs(input,map)
      % Contruct grid.
      default = struct('type',string.empty);
      params = readInput(default,input);

      data = zeros(map);
      switch lower(params.type)
        case "all"
          data(:)=1:prod(map);
        case "surface"
          % TO DO implement this input
        otherwise
          gresLog().error("The grid parameters for the simulation is not well defined!");
      end
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






  end

end