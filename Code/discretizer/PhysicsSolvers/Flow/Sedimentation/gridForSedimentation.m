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

    function [id,faceArea, dh, matP] = getBordCell(obj,label)
      % GETBORDCELL Returns IDs, areas, and half-heights for a boundary.
      % Labels: 'x0', 'xm', 'y0', 'ym', 'z0' (bottom), 'zm' (top).

      nmat = size(obj.matfrac,2);
      switch lower(label)
        case "x0"
          areas = diff(obj.grid.Y).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,1);
          id=zeros(nelm(1),1,"uint64");
          faceArea=zeros(nelm(1),1);
          dh=zeros(nelm(1),1);
          count=1;
          for j=1:obj.ncells(2)
            for h=1:obj.columnsHeight(1,j)
              id(count)=obj.dof(1,j,h);
              faceArea(count)=areas(j,h);
              count=count+1;
            end
          end
        case "xm"
          areas = diff(obj.grid.Y).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,1);
          id=zeros(nelm(end),1,"uint64");
          faceArea=zeros(nelm(end),1);
          dh=zeros(nelm(1),1);
          count=1;
          for j=1:obj.ncells(2)
            for h=1:obj.columnsHeight(end,j)
              id(count)=obj.dof(end,j,h);
              faceArea(count)=areas(j,h);
              count=count+1;
            end
          end
        case "y0"
          areas = diff(obj.grid.X).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,2);
          id=zeros(nelm(1),1,"uint64");
          faceArea=zeros(nelm(1),1);
          dh=zeros(nelm(1),1);
          count=1;
          for i=1:obj.ncells(1)
            for h=1:obj.columnsHeight(i,1)
              id(count)=obj.dof(i,1,h);
              faceArea(count)=areas(i,h);
              count=count+1;
            end
          end
        case "ym"
          areas = diff(obj.grid.X).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,2);
          id=zeros(nelm(end),1,"uint64");
          faceArea=zeros(nelm(end),1);
          dh=zeros(nelm(1),1);
          count=1;
          for i=1:obj.ncells(1)
            for h=1:obj.columnsHeight(i,end)
              id(count)=obj.dof(i,end,h);
              faceArea(count)=areas(i,h);
              count=count+1;
            end
          end
        case "z0"
          areas = diff(obj.grid.Y).*diff(obj.grid.X)';
          nelm = prod(obj.ncells(1:2));
          id=zeros(nelm,1,"uint64");
          faceArea=zeros(nelm,1);
          dh(1:nelm,1)=(obj.grid.Z(2)-obj.grid.Z(1))/2;
          count=1;
          for i=1:obj.ncells(1)
            for j=1:obj.ncells(2)
              id(count)=obj.dof(i,j,1);
              faceArea(count)=areas(i,j);
              count=count+1;
            end
          end
        case "zm"
          dhtmp=diff(obj.grid.Z)/2;
          areas = diff(obj.grid.Y).*diff(obj.grid.X)';
          nelm = prod(obj.ncells(1:2));
          id=zeros(nelm,1,"uint64");
          faceArea=zeros(nelm,1);
          dh=zeros(nelm,1);
          count=1;
          for i=1:obj.ncells(1)
            for j=1:obj.ncells(2)
              dh(count)=dhtmp(obj.columnsHeight(i,j));
              id(count)=obj.dof(i,j,obj.columnsHeight(i,j));
              faceArea(count)=areas(i,j);
              count=count+1;
            end
          end
        otherwise
          id = [];
          faceArea = [];
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
      ncell = double(obj.ncells);
      cellID = uint64(prod(ncell(1:2))*(k-1)+ncell(1)*(j-1)+i);
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


    function grow(obj,cells)
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
      obj.columnsHeight = reshape(obj.columnsHeight,obj.ncells(1:2));

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



  end

end