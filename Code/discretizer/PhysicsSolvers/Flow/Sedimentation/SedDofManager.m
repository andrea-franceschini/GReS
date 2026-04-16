classdef SedDofManager < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties
    % Dynamic grid control
    ncells (1,3)                 % Number of cells in (x,y,z)
    columnsHeight (:,:)          % Active height per column
    dof (:,:,:)                  % Cell DOF map (i,j,k) -> DOF ID
    ndofs                        % Number of active DOFs
  end

  methods (Access = public)
    function obj = SedDofManager(data)
      obj.ncells = data.Grid.division;
      obj.dof = zeros(obj.ncells);
      obj.columnsHeight = zeros(prod(obj.ncells(1:2)),1);

      % Constructing the initial grid
      switch data.Initial.type
        case "all"
          obj.dof(:) = 1:prod(obj.ncells);
          obj.columnsHeight(:) = obj.ncells(3);
        case "surface"
          numcells = prod(obj.ncells(1:2));

          if isfield(data,"columnHeight")
            surfCellHeight = load(data.columnHeight);
          end

          count=1;
          for j=1:obj.ncells(2)
            for i=1:obj.ncells(1)
              obj.dof(i,j,surfCellHeight(i,j))=count;
              count=count+1;
            end
          end

          obj.columnsHeight = surfCellHeight(:);
        case "untilSurface"
        case "doubleSurface"
        otherwise
          gridForSedimentation.errorGrid();
      end

      obj.ndofs = sum(obj.dof~=0,"all");
    end

    function dofs = getBordXMin(obj)
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
      dofs = getBordX(obj,J_idx,nelm);
    end

    function dofs = getBordXMax(obj)
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
      dofs = getBordX(obj,J_idx,nelm);
    end

    function dofs = getBordYMin(obj)
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
      dofs = getBordY(obj,I_idx,nelm);
    end

    function dofs = getBordYMax(obj)
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
      dofs = getBordY(obj,I_idx,nelm);
    end

    function dofs = getBordZMin(obj)
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = ones(prod(obj.ncells(1:2)),1);

      % For case where some column do not exist
      % map = obj.columnsHeight~=0;
      % cellID = sub2ind(obj.ncells,idI(map),idJ(map),idK(map));

      cellID = sub2ind(obj.ncells,idI,idJ,idK);
      dofs = obj.dof(cellID);
    end

    function dofs = getBordZMax(obj)
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = obj.columnsHeight;

      % For case where some column do not exist
      % map = obj.columnsHeight~=0;
      % cellID = sub2ind(obj.ncells,idI(map),idJ(map),idK(map));

      cellID = sub2ind(obj.ncells,idI,idJ,idK);
      dofs = obj.dof(cellID);
    end

    function neigh = getNeigh(obj,dofs)
      % GETNEIGH Returns neighboring cell DOFs.
      %
      % Order:
      %   [x-, x+, y-, y+, z-, z+]

      if ~exist("dofs","var")
        dofs = (1:obj.ndofs)';
      end

      % pos = find(ismember(obj.dof,dofs));
      map = ismember(obj.dof,dofs);
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [idI,idJ,idK]=ind2sub(obj.ncells,pos(id));
      neigh = zeros(length(dofs),6);

      actTmp = idJ~=1;
      pos = sub2ind(obj.ncells, idI(actTmp), idJ(actTmp)-1, idK(actTmp));
      neigh(actTmp,1) =  obj.dof(pos);
      actTmp = idJ~=obj.ncells(2);
      pos = sub2ind(obj.ncells, idI(actTmp), idJ(actTmp)+1, idK(actTmp));
      neigh(actTmp,2) =  obj.dof(pos);

      actTmp = idI~=1;
      pos = sub2ind(obj.ncells, idI(actTmp)-1, idJ(actTmp), idK(actTmp));
      neigh(actTmp,3) =  obj.dof(pos);
      actTmp = idI~=obj.ncells(1);
      pos = sub2ind(obj.ncells, idI(actTmp)+1, idJ(actTmp), idK(actTmp));
      neigh(actTmp,4) =  obj.dof(pos);      

      actTmp = idK~=1;
      pos = sub2ind(obj.ncells, idI(actTmp), idJ(actTmp), idK(actTmp)-1);
      neigh(actTmp,5) =  obj.dof(pos);
      actTmp = idK~=obj.ncells(3);
      pos = sub2ind(obj.ncells, idI(actTmp), idJ(actTmp), idK(actTmp)+1);
      neigh(actTmp,6) =  obj.dof(pos);
    end

    function out = getMapFromDofs(obj,dofs)  %<--- ok
      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      map = ismember(obj.dof,dofs);
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [idI,idJ,~]=ind2sub(obj.ncells,pos(id));
      out = sub2ind(obj.ncells(1:2), idI', idJ');
    end

    function dof = getMaxDofUnchanged(obj)
      minval = min(obj.columnsHeight);
      if minval~=1
        loc = obj.columnsHeight == minval;
        [idI, idJ] = obj.getIJLay();
        pos = sub2ind(obj.ncells, idI(loc), idJ(loc), obj.columnsHeight(loc)-1);
        dof = max(obj.dof(pos));
      else
        dof = 0;
      end
    end

    function [idI,idJ,idK] = getIJKFromDofs(obj,dofs)
      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      map = ismember(obj.dof,dofs);
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [idI,idJ,idK]=ind2sub(obj.ncells,pos(id));
    end

    function newlayer = grow(obj,map)
      % GROW Adds new cells on top of selected columns.
      %
      % Inputs:
      %   map      - Logical array of growing columns
      %   height   - Height of the new layer
      %
      % Output:
      %   newlayer - True if a new Z layer was added

      cellsTadd = sum(map);

      % check if is necessary to add a new layer in z
      atTop = obj.columnsHeight == obj.ncells(3);
      newlayer = any(and(atTop,map));
      if newlayer
        obj.ncells(3) = obj.ncells(3)+1;
        obj.dof(:,:,end+1) = zeros(obj.ncells(1:2));
      end

      % update the dof for the newest cells.
      [idI,idJ] = obj.getIJLay;
      pos = sub2ind(obj.ncells,idI(map),idJ(map),obj.columnsHeight(map)+1);
      obj.dof(pos) = (obj.ndofs+1:obj.ndofs+cellsTadd)';
      obj.ndofs = obj.ndofs + cellsTadd;
      obj.columnsHeight(map) = obj.columnsHeight(map)+1;
    end

    function npts = getNumberPoints(obj)
      % GETNUMBERPOINTS Returns total number of grid points.
      npts = (obj.ncells(1)+1)*(obj.ncells(2)+1)*(obj.ncells(3)+1);
    end


    function data = accByColumnFromBot2Top(obj,valByCell)
      data = zeros(obj.ndofs,1);
      for i=1:obj.ncells(1)
        for j=1:obj.ncells(2)
          column = obj.dof(i,j,:);
          acc = 0.;
          for k=1:obj.ncells(3)
            dofId = column(k);
            if dofId ~= 0
              acc = acc + valByCell(dofId);
              data(dofId)=acc;
            end
          end
        end
      end
    end

    function data = cell2NodeAccByColumnFromBot2Top(obj,valByCell)
      data = zeros(8*obj.ndofs,1);
      for i=1:obj.ncells(1)
        for j=1:obj.ncells(2)
          column = obj.dof(i,j,:);
          acc = 0.;
          for k=1:obj.ncells(3)
            dofId = column(k);
            if dofId ~= 0
              refDof = 8*(dofId-1);
              data(refDof+1:refDof+4)=acc;
              acc = acc + valByCell(dofId);
              data(refDof+5:refDof+8)=acc;              
            end
          end
        end
      end
    end


  end


  methods (Access = private)

    function id = getBordX(obj,idx,nelm)
      id=zeros(nelm,1);

      count=1;
      for i=1:obj.ncells(1)
        for k=1:obj.ncells(3)
          pos = obj.ncells(1)*(i-1)+idx(i);
          if k>obj.columnsHeight(pos)
            break;
          end
          id(count)=obj.dof(i,idx(i),k);
          count=count+1;
        end
      end
    end

    function id = getBordY(obj,idx,nelm)
      id=zeros(nelm,1);

      count=1;
      for j=1:obj.ncells(2)
        for k=1:obj.ncells(3)
          pos = obj.ncells(1)*(idx(j)-1)+j;
          if k>obj.columnsHeight(pos)
            break;
          end
          id(count)=obj.dof(idx(j),j,k);
          count=count+1;
        end
      end
    end

    function [idI, idJ] = getIJLay(obj)
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
    end



  end
end