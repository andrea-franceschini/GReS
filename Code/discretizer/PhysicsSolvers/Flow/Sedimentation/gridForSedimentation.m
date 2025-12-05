classdef gridForSedimentation < handle
  %SEDGRID Class to control the grid in the sedimentation physics
  %   Detailed explanation goes here

  properties
    grid = structGrid()
    ncells (1,3) int64
    columnsHeight (:,:) int64
    mapCellIds (:,:,:) int64
    mapMatFlags (:,:,:) int64
  end

  methods (Access = public)
    function obj = gridForSedimentation(varargin)

      if nargin == 0
        return
      end
      data = varargin{1};
      obj.grid = structGrid(data.grid);
      obj.ncells = obj.grid.getNumberCells;

      obj.columnsHeight = zeros(obj.ncells(1:2));
      obj.mapCellIds = zeros(obj.ncells);
      obj.mapMatFlags = zeros(obj.ncells);
      for i=1:length(data.initial)
        obj.mapMatFlags(:,:,i)=data.initial(i).lay.materialFlag;
        tmp=max(obj.mapCellIds,[],"all");
        obj.mapCellIds(:,:,i)=reshape(tmp+1:prod(obj.ncells(1:2)),obj.ncells(1:2));
      end
      obj.columnsHeight(:) = length(data.initial);

    end

    function id = getBordCellID(obj,label)
      % Function the get the number of the cell in the border.
      switch lower(label)
        case "x0"
          id = obj.mapCellIds(1,:,:);
        case "xm"
          id = obj.mapCellIds(end,:,:);
        case "y0"
          id = obj.mapCellIds(:,1,:);
        case "ym"
          id = obj.mapCellIds(:,end,:);
        case "z0"
          id = obj.mapCellIds(:,:,1);
        case "zm"
          id = zeros(obj.ncells(1:2));
          for j=1:obj.ncells(2)
            for i=1:obj.ncells(1)
              id(i,j)=obj.columnsHeight(i,j);
            end
          end
        otherwise
          id = [];
          return
      end
      id = id(:);
    end


  end
end