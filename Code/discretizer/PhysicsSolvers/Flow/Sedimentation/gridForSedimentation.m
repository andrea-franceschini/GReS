classdef gridForSedimentation < handle
  %SEDGRID Class to control the grid in the sedimentation physics
  %   Detailed explanation goes here

  properties
    % Grid information
    grid = structGrid()
    ncells (1,3) uint64

    % Control for the material
    matfrac (:,:) % Fraction of each material per cell (cells,nmat)

    % Control for the growing grid
    columnsHeight (:,:) uint64
    mapCellIds (:,:,:) uint64
    numberActiveCells uint64
    


  end

  properties(Access = private)

  end

  methods (Access = public)
    function obj = gridForSedimentation(varargin)
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
      if flagEscape
        return
      end

      % Assuming only 1 material
      if ~exist('matLabel', 'var') == 1
        matLabel = 1;
      end

      obj.constructor(data,matLabel);
    end

    function [id,faceArea, dh, matP] = getBordCell(obj,label)
      % Function the get the number of the cell in the border.

      nmat = size(obj.matfrac,2);
      switch lower(label)
        case "x0"
          areas = diff(obj.grid.Y).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,1);
          id=zeros(nelm(1),1,"uint64");
          faceArea=zeros(nelm(1),1);
          % matP=zeros(nelm(1),nmat);
          dh=zeros(nelm(1),1);
          count=1;
          for j=1:obj.ncells(2)
            for h=1:obj.columnsHeight(1,j)
              id(count)=obj.mapCellIds(1,j,h);
              faceArea(count)=areas(j,h);
              % matP(count,:)=obj.matfrac(id(count),:);
              count=count+1;
            end
          end
          % id = obj.mapCellIds(1,:,:);
        case "xm"
          areas = diff(obj.grid.Y).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,1);
          id=zeros(nelm(end),1,"uint64");
          faceArea=zeros(nelm(end),1);
          dh=zeros(nelm(1),1);
          count=1;
          for j=1:obj.ncells(2)
            for h=1:obj.columnsHeight(end,j)
              id(count)=obj.mapCellIds(end,j,h);
              faceArea(count)=areas(j,h);
              count=count+1;
            end
          end
          % id = obj.mapCellIds(end,:,:);
        case "y0"
          areas = diff(obj.grid.X).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,2);
          id=zeros(nelm(1),1,"uint64");
          faceArea=zeros(nelm(1),1);
          dh=zeros(nelm(1),1);
          count=1;
          for i=1:obj.ncells(1)
            for h=1:obj.columnsHeight(i,1)
              id(count)=obj.mapCellIds(i,1,h);
              faceArea(count)=areas(i,h);
              count=count+1;
            end
          end
          % id = obj.mapCellIds(:,1,:);
        case "ym"
          areas = diff(obj.grid.X).*diff(obj.grid.Z)';
          nelm = sum(obj.columnsHeight,2);
          id=zeros(nelm(end),1,"uint64");
          faceArea=zeros(nelm(end),1);
          dh=zeros(nelm(1),1);
          count=1;
          for i=1:obj.ncells(1)
            for h=1:obj.columnsHeight(i,end)
              id(count)=obj.mapCellIds(i,end,h);
              faceArea(count)=areas(i,h);
              count=count+1;
            end
          end
          % id = obj.mapCellIds(:,end,:);
        case "z0"
          areas = diff(obj.grid.Y).*diff(obj.grid.X)';
          nelm = prod(obj.ncells(1:2));
          id=zeros(nelm,1,"uint64");
          faceArea=zeros(nelm,1);
          dh(1:nelm,1)=(obj.grid.Z(2)-obj.grid.Z(1))/2;
          count=1;
          for i=1:obj.ncells(1)
            for j=1:obj.ncells(2)
              id(count)=obj.mapCellIds(i,j,1);
              faceArea(count)=areas(i,j);
              count=count+1;
            end
          end
          % id = obj.mapCellIds(:,:,1);
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
              id(count)=obj.mapCellIds(i,j,obj.columnsHeight(i,j));
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
      ncells = sum(sum(obj.columnsHeight));
    end




    function coord = getCoordCenter(obj,cellIds)
      % Function to return the size of the block

      % Return the coordinates from the layer kpos
      tmp = zeros(obj.numberActiveCells,3);
      count=1;
      for k=1:obj.ncells(3)
        ztmp = obj.grid.centerZ(k);
        for j=1:obj.ncells(2)
          ytmp = obj.grid.centerY(j);
          for i=1:obj.ncells(1)
            tmp(count,1)=obj.grid.centerX(i);
            tmp(count,2)=ytmp;
            tmp(count,3)=ztmp;
            count=count+1;
          end
        end
      end
      coord=tmp(cellIds,:);
    end



    function vols = computeVols(obj)
      segmX = reshape(diff(obj.grid.X),3,1,1);
      segmY = reshape(diff(obj.grid.Y),1,3,1);
      segmZ = reshape(diff(obj.grid.Z),1,1,3);
      vols=segmX .* segmY .* segmZ;
      cellsOn = obj.mapCellIds ~= 0;
      map = obj.mapCellIds(cellsOn);
      vols=vols(cellsOn);
      vols=vols(map);
    end


    function [i,j,k] = getIJKfromCellID(obj,cellID)
      [i,j,k]=ind2sub(obj.ncells,cellID);
    end

    function cellID = getCellIDfromIJK(obj,i,j,k)
      ncell = double(obj.ncells);
      surf=prod(ncell(1:2));
      cellID = uint64(surf*(k-1)+ncell(1)*(j-1)+i);
    end

    function cellID = getActiveCells(obj)
      cellID = obj.mapCellIds(obj.mapCellIds ~=0);
    end
   


  end

  methods (Access = private)
    function constructor(obj,data,matList)
      obj.grid = structGrid(data.grid);
      obj.ncells = obj.grid.getNumberCells;

      nmat = length(matList);
      obj.matfrac = zeros(prod(obj.ncells),nmat);

      obj.columnsHeight = zeros(obj.ncells(1:2));
      obj.mapCellIds = zeros(obj.ncells);

      cellsByLay = (1:prod(obj.ncells(1:2)))';
      for i=1:length(data.initial.lay)
        lay = data.initial.lay(i);

        % Cells referring to this layer
        if isfield(lay,"number")
          posInitial = prod(obj.ncells(1:2))*(lay.number-1);
        end
        if isfield(lay,"materialFlag")
          mats = double(matList == lay.materialFlag)';
        end
        if isfield(lay,"fractions")
          mats = double(lay.fractions)';
        end
        obj.matfrac(cellsByLay+posInitial,:) = mats.*ones(prod(obj.ncells(1:2)),nmat);

        tmp=max(obj.mapCellIds,[],"all");
        obj.mapCellIds(:,:,i)=reshape(tmp+1:tmp+prod(obj.ncells(1:2)),obj.ncells(1:2));
      end
      obj.columnsHeight(:) = length(data.initial.lay);

      obj.numberActiveCells = sum(obj.columnsHeight,"all");
    end


    
  end

  methods (Static)
    
  end
end