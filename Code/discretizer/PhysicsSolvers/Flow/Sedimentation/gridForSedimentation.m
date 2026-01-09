classdef gridForSedimentation < handle
  % GRIDFORSEDIMENTATION Manages a dynamic 3D grid for sedimentation.
  %   Tracks material fractions, active cells, and variable column heights
  %   to simulate growth by material deposition.

  properties
    % Grid information
    grid = structGrid()
    ncells (1,3) uint64

    % Control for the material
    matfrac (:,:) % Fraction of each material per cell (cells,nmat)

    % Control for the growing grid
    columnsHeight (:,:) uint64
    dof (:,:,:) uint64
    ndofs uint64
  end

  properties(Access = private)
    nmat uint16
  end

  methods (Access = public)
    function obj = gridForSedimentation(varargin)
      % Constructor: Initialize via XML data and material labels.
      flagEscape=true;
      for k = 1:2:nargin
        switch lower(varargin{k})
          case "xml"
            data = varargin{k+1};
            flagEscape = false;
          case "materiais"
            matLabel = varargin{k+1};
          case "nummateriais"
            matLabel = (1:varargin{k+1})';
        end
      end

      % Check if the class is gonna be constructed.
      if flagEscape, return; end

      % Assuming only 1 material
      if ~exist('matLabel', 'var'), matLabel = 1; end

      obj.constructor(data,matLabel);
    end

    function [id,faceArea,dh] = getBordCell(obj,label)
      % GETBORDCELL Returns IDs, areas, and half-heights for a boundary.
      % Labels: 'x0', 'xm', 'y0', 'ym', 'z0' (bottom), 'zm' (top).

      switch lower(label)
        case "x0"
          mapH = reshape(obj.columnsHeight,obj.ncells(1:2));
          cells = mapH(:,1);
          J_idx = zeros(obj.ncells(1),1);
          for i=1:obj.ncells(1)
            count=1;
            while mapH(i,count) == 0 && count<=obj.ncells(2)
              count=count+1;
            end
            cells(i)=mapH(i,count);
            J_idx(i)=count;
          end
          nelm = sum(cells,"all");
          [id,faceArea,dh] = getBordX(obj,J_idx,nelm);

        case "xm"
          mapH = reshape(obj.columnsHeight,obj.ncells(1:2));
          cells = mapH(:,end);
          J_idx = zeros(obj.ncells(1),1);
          for i=1:obj.ncells(1)
            count=obj.ncells(2);
            while mapH(i,count) == 0 && count>=1
              count=count-1;
            end
            cells(i)=mapH(i,count);
            J_idx(i)=count;
          end
          nelm = sum(cells,"all");
          [id,faceArea,dh] = getBordX(obj,J_idx,nelm);

        case "y0"
          mapH = reshape(obj.columnsHeight,obj.ncells(1:2));
          cells = mapH(1,:);
          I_idx = zeros(obj.ncells(2),1);
          for j=1:obj.ncells(2)
            count=1;
            while mapH(count,j) == 0 && count<=obj.ncells(2)
              count=count+1;
            end
            cells(j)=mapH(count,j);
            I_idx(j)=count;
          end
          nelm = sum(cells,"all");
          [id,faceArea,dh] = getBordY(obj,I_idx,nelm);

        case "ym"
          mapH = reshape(obj.columnsHeight,obj.ncells(1:2));
          cells = mapH(end,:);
          I_idx = zeros(obj.ncells(2),1);
          for j=1:obj.ncells(2)
            count=obj.ncells(1);
            while mapH(count,j) == 0 && count>=1
              count=count-1;
            end
            cells(j)=mapH(count,j);
            I_idx(j)=count;
          end
          nelm = sum(cells,"all");
          [id,faceArea,dh] = getBordY(obj,I_idx,nelm);

        case "z0"
          [id,faceArea,dh] = getBordZ(obj,true);

        case "zm"
          [id,faceArea,dh] = getBordZ(obj,false);

        otherwise
          id = [];
          faceArea = [];
          dh = [];
          return
      end
    end

    function ncells = getNumberCells(obj)
      % Returns total active cell count.
      ncells = sum(obj.columnsHeight, 'all');
    end

    function coord = getCoordCenter(obj, cellIds)
      % GETCOORDCENTER Returns XYZ center coordinates for specified cell IDs.
      %   Uses ndgrid for vectorized coordinate generation to improve speed.

      % Generate 3D matrices for each coordinate component
      [X, Y, Z] = ndgrid(obj.grid.centerX, obj.grid.centerY, obj.grid.centerZ);

      % Extract the coordinates for the requested linear indices
      % This assumes cellIds corresponds to standard (i,j,k) linear indexing
      activeIds = obj.getActiveCells();
      coord = [X(activeIds), Y(activeIds), Z(activeIds)];
      coord = coord(obj.dof(activeIds),:);
      coord = coord(cellIds,:);
    end

    function vols = computeVols(obj)
      % COMPUTEVOLS Calculates volumes for all active cells using broadcasting.
      dx = diff(obj.grid.X(:));
      dy = diff(obj.grid.Y(:));
      dz = diff(obj.grid.Z(:));

      % Create 3D volume matrix via outer product (broadcasting)
      vAll = dx .* reshape(dy, 1, []) .* reshape(dz, 1, 1, []);

      % Extract volumes only for active cells
      activeIds = obj.getActiveCells();
      vols = vAll(activeIds);
      vols = vols(obj.dof(activeIds));
    end

    function ijk = getIJKfromCellID(obj,cellID)
      % Converts linear Cell ID to 3D subscripts.
      [i,j,k]=ind2sub(obj.ncells,cellID);
      ijk=[i,j,k];
    end

    function cellID = getCellIDfromIJK(obj,i,j,k)
      % Converts 3D subscripts to linear Cell ID.
      cellID = sub2ind(obj.ncells,i,j,k);
      % ncell = double(obj.ncells);
      % cellID = uint64(prod(ncell(1:2))*(k-1)+ncell(1)*(j-1)+i);
    end

    function [idI,idJ,idK] = getIJKTop(obj)
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = obj.columnsHeight;
    end

    function cellID = getActiveCells(obj)
      % Returns a list of all non-zero Cell IDs.
      cellID = obj.dof(obj.dof ~=0);
    end

    function dofs = getActiveDofs(obj)
      dofs = obj.dof(obj.dof ~=0);
    end

    function ndofs = getNdofs(obj)
      ndofs=obj.ndofs;
    end

    function dofs = getDofsFromIJK(obj,ijk)
      if nargin == 1
        dofs = obj.dof(:);
      else
        pos = sub2ind(size(obj.dof), ijk(:,1), ijk(:,2), ijk(:,3));
        dofs = obj.dof(pos);
      end
    end

    function neigh = getNeigh(obj,ijk)
      neigh = zeros(size(ijk,1),6);

      actTmp = ijk(:,1)~=1;
      pos = sub2ind(size(obj.dof), ijk(actTmp,1)-1, ijk(actTmp,2), ijk(actTmp,3));
      neigh(actTmp,3) =  obj.dof(pos);
      actTmp = ijk(:,1)~=obj.ncells(1);
      pos = sub2ind(size(obj.dof), ijk(actTmp,1)+1, ijk(actTmp,2), ijk(actTmp,3));
      neigh(actTmp,4) =  obj.dof(pos);

      actTmp = ijk(:,2)~=1;
      pos = sub2ind(size(obj.dof), ijk(actTmp,1), ijk(actTmp,2)-1, ijk(actTmp,3));
      neigh(actTmp,1) =  obj.dof(pos);
      actTmp = ijk(:,2)~=obj.ncells(2);
      pos = sub2ind(size(obj.dof), ijk(actTmp,1), ijk(actTmp,2)+1, ijk(actTmp,3));
      neigh(actTmp,2) =  obj.dof(pos);

      actTmp = ijk(:,3)~=1;
      pos = sub2ind(size(obj.dof), ijk(actTmp,1), ijk(actTmp,2), ijk(actTmp,3)-1);
      neigh(actTmp,5) =  obj.dof(pos);
      actTmp = ijk(:,3)~=obj.ncells(3);
      pos = sub2ind(size(obj.dof), ijk(actTmp,1), ijk(actTmp,2), ijk(actTmp,3)+1);
      neigh(actTmp,6) =  obj.dof(pos);
    end

    function newlayer = grow(obj,map,matfrac,height)
      cellsTadd = size(matfrac,1);    

      % check if is necessary to add a new layer in z
      atTop = obj.columnsHeight == obj.ncells(3);
      newlayer = any(and(atTop,map));
      if newlayer
        % update the structured grid
        addCoordZ(obj.grid,height);
        obj.dof(:,:,end+1) = zeros(obj.ncells(1:2),"uint64");
        obj.ncells(3)=obj.ncells(3)+1;
      end

      obj.matfrac(end+1:end+cellsTadd,:) = matfrac;
      pos = sub2ind(obj.ncells, ...
        repmat((1:obj.ncells(1))',obj.ncells(2),1), ...
        repelem((1:obj.ncells(2))',obj.ncells(1)), ...
        obj.columnsHeight+1);
      pos = pos(map);
      obj.dof(pos) = (obj.ndofs+1:obj.ndofs+cellsTadd)';
      obj.ndofs = obj.ndofs + cellsTadd;
      obj.columnsHeight(map)=obj.columnsHeight(map)+1;
    end


  end

  methods (Access = private)
    function constructor(obj,data,matList)
      % Internal initialization of grid, maps, and material layers.
      obj.grid = structGrid(data.grid);
      obj.ncells = obj.grid.getNumberCells;
      obj.ndofs = prod(obj.ncells);

      nCellMap = prod(obj.ncells(1:2));

      obj.nmat = length(matList);
      obj.matfrac = zeros(obj.ndofs,obj.nmat);
      obj.dof = reshape(1:obj.ndofs,obj.ncells);
      obj.columnsHeight = obj.ncells(3)*ones(nCellMap,1,"uint64");

      for lay=data.initial.all
        if isfield(lay,"materialFlag")
          obj.matfrac(:,lay.materialFlag)=1.;
        end
      end

      cellsByLay = (1:nCellMap)';
      for lay=data.initial.lay
        % Cells referring to this layer
        if isfield(lay,"number")
          posInitial = nCellMap*(lay.number-1);
          pos = cellsByLay+posInitial;
        end
        if isfield(lay,"materialFlag")
          mats = double(matList == lay.materialFlag)';
        end
        if isfield(lay,"fractions")
          mats = double(lay.fractions)';
        end
        obj.matfrac(pos,:)=mats.*ones(nCellMap,obj.nmat);
      end
    end


    function [id,area,dh] = getBordX(obj,idx,nelm)
      dhtmp=diff(obj.grid.Z)/2;
      areas = diff(obj.grid.X).*diff(obj.grid.Z)';

      id=zeros(nelm,1,"uint64");
      area=zeros(nelm,1);
      dh=zeros(nelm,1);

      count=1;
      for i=1:obj.ncells(1)
        for k=1:obj.ncells(3)
          pos = obj.ncells(1)*(i-1)+idx(i);
          if k>obj.columnsHeight(pos)
            break;
          end
          dh(count)=dhtmp(k);
          id(count)=obj.dof(i,idx(i),k);
          area(count)=areas(i,k);
          count=count+1;
        end
      end
    end


    function [id,area,dh] = getBordY(obj,idx,nelm)
      dhtmp=diff(obj.grid.Z)/2;
      areas = diff(obj.grid.Y).*diff(obj.grid.Z)';

      id=zeros(nelm,1,"uint64");
      area=zeros(nelm,1);
      dh=zeros(nelm,1);

      count=1;
      for j=1:obj.ncells(2)
        for k=1:obj.ncells(3)
          pos = obj.ncells(1)*(idx(j)-1)+j;
          if k>obj.columnsHeight(pos)
            break;
          end
          dh(count)=dhtmp(k);
          id(count)=obj.dof(idx(j),j,k);
          area(count)=areas(idx(j),k);
          count=count+1;
        end
      end
    end


    function [id,area,dh] = getBordZ(obj,bord)
      activeCells = obj.columnsHeight>0 ;
      nelm = sum(activeCells);

      dhtmp=diff(obj.grid.Z)/2;
      areas = reshape(diff(obj.grid.Y).*diff(obj.grid.X)',nelm,1);

      id=zeros(nelm,1,"uint64");
      area=zeros(nelm,1);
      dh=zeros(nelm,1);

      count=1;
      if bord
        % Case for the bord at z0
        for column=1:length(activeCells)
          if activeCells(column)
            I_idx = mod(column-1,obj.ncells(1))+1;
            J_idx = (column-I_idx)/obj.ncells(1)+1;
            dh(count)=dhtmp(1);
            id(count)=obj.dof(I_idx,J_idx,1);
            area(count)=areas(column);
          end
          count=count+1;
        end
      else
        % Case for the bord at zm
        for column=1:length(activeCells)
          if activeCells(column)
            I_idx = mod(column-1,obj.ncells(1))+1;
            J_idx = (column-I_idx)/obj.ncells(1)+1;
            dh(count)=dhtmp(obj.columnsHeight(column));
            id(count)=obj.dof(I_idx,J_idx,obj.columnsHeight(column));
            area(count)=areas(column);
          end
          count=count+1;
        end
      end
    end



  end

end