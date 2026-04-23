classdef SedDofManager < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties
    % Dynamic grid control
    ncells (1,3)                 % Number of cells in (x,y,z)
    laysByCol (:,:)              % Number of layer by column
    dof (:,:,:)                  % Cell DOF map (i,j,k) -> DOF ID
    ndofs                        % Number of active DOFs
    inactDof                     % inactive dof
  end

  properties (Access = private)
    precs = 1e-5
  end

  methods (Access = public)
    function obj = SedDofManager(data,ncells)
      obj.ncells = ncells;
      obj.dof = zeros(obj.ncells);
      % obj.laysByCol = zeros(prod(obj.ncells(1:2)),1);
      obj.laysByCol = zeros(obj.ncells(1:2));

      % Constructing the initial grid
      switch data.Initial.type
        case "all"
          obj.dof(:) = 1:prod(obj.ncells);
          obj.laysByCol(:) = obj.ncells(3);
          obj.inactDof = (prod(ncells(1:2))+1:prod(ncells))';
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

          obj.laysByCol = surfCellHeight(:);
        case "untilSurface"
        case "doubleSurface"
        otherwise
          gridForSedimentation.errorGrid();
      end

      obj.ndofs = sum(obj.dof~=0,"all");
    end

    function dofs = getActiveDof(obj)
      dofs = 1:obj.ndofs;
      dofs = dofs(~ismember(dofs,obj.inactDof));
    end

    function dofs = getNotActiveDof(obj)
      dofs = obj.inactDof;
    end

    function dofs = getActiveNdof(obj)
      dofs = obj.ndofs-length(obj.inactDof);
    end

    function dofs = activetedDofs(obj,map)
      % Activeted the inactive dofs using the map as reference.
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = obj.findIdLaysVaryCell;

      cellID = sub2ind(obj.ncells,idI(map(:)),idJ(map(:)),idK(map(:)));
      actDof = obj.dof(cellID);
      loc = ismember(obj.inactDof,actDof);
      dofs = obj.inactDof(loc);
      obj.inactDof = obj.inactDof(~loc);
    end

    function dofs = getVarHeightDofs(obj,map)
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = obj.findIdLaysVaryCell;
      cellID = sub2ind(obj.ncells,idI(map),idJ(map),idK(map(:)));
      dofs = obj.dof(cellID);
    end

    function deactiveDofs(obj,colDeact)
      dofs = obj.getVarHeightDofs(colDeact);
      ncell = length(dofs);
      obj.inactDof(end+1:end+ncell) = dofs;
    end


    function mesh = makeMeshOutput(obj,mesh,coordX,coordY,coordZ,sedAcc,act)
      idxTop = obj.findIdLaysVaryCell;
      idxTop(~act) = idxTop(~act)-1;

      % Remap the dofs
      allDofs = 1:obj.ndofs;
      actDofs = obj.getActiveDof;
      [tf, loc] = ismember(allDofs, actDofs);
      dofs = loc;
      dofs(~tf) = 0;

      npoint = sum(4*idxTop+4,"all");
      ncell = sum(idxTop,"all");

      coord(1:npoint,1:3) = 0.;
      conec(1:ncell,1:8)  = 0;

      count = 0;
      for col=1:prod(obj.ncells(1:2))
        idI = mod(col-1,obj.ncells(1))+1;
        idJ = (col-idI)/obj.ncells(1)+1;

        coord(4*count+1,1) = coordX(idI);
        coord(4*count+2,1) = coordX(idI+1);
        coord(4*count+3,1) = coordX(idI+1);
        coord(4*count+4,1) = coordX(idI);
        
        coord(4*count+1,2) = coordY(idJ);
        coord(4*count+2,2) = coordY(idJ);
        coord(4*count+3,2) = coordY(idJ+1);
        coord(4*count+4,2) = coordY(idJ+1);
        
        coord(4*count+1:4*count+4,3) = coordZ(1);
        count = count+1;
        topLay = idxTop(col);
        for lay = 1:topLay
          coord(4*count+1,1) = coordX(idI);
          coord(4*count+2,1) = coordX(idI+1);
          coord(4*count+3,1) = coordX(idI+1);
          coord(4*count+4,1) = coordX(idI);

          coord(4*count+1,2) = coordY(idJ);
          coord(4*count+2,2) = coordY(idJ);
          coord(4*count+3,2) = coordY(idJ+1);
          coord(4*count+4,2) = coordY(idJ+1);

          if lay == topLay
            if act(col)
              coord(4*count+1:4*count+4,3) = coordZ(lay)+sedAcc(col);
            else
              coord(4*count+1:4*count+4,3) = coordZ(lay+1);
            end
          else
            coord(4*count+1:4*count+4,3) = coordZ(lay+1);
          end
          cell = dofs(obj.dof(idI,idJ,lay));
          conec(cell,1)=4*(count-1)+1;
          conec(cell,2)=4*(count-1)+2;
          conec(cell,3)=4*(count-1)+3;
          conec(cell,4)=4*(count-1)+4;
          conec(cell,5)=4*count+1;
          conec(cell,6)=4*count+2;
          conec(cell,7)=4*count+3;
          conec(cell,8)=4*count+4;
          count=count+1;
        end
      end

      mesh.nDim = 3;
      mesh.nCells = ncell;      
      mesh.cellVTKType = 12*ones(ncell,1);
      mesh.cellNumVerts = 8*ones(ncell,1);
      mesh.cellTag = 8*ones(ncell,1);

      mesh.nNodes = npoint;
      mesh.coordinates = coord;
      mesh.cells = conec;
    end

    function data = getComp(obj,act,comp)
      % Identify the top layer index.
      % sedAcc = reshape(sedAcc,obj.ncells(1:2));
      % sedMin = control*obj.precs;
      % act = sedAcc > sedMin;
      idxTop = obj.findIdLaysVaryCell;
      idxTop(~act) = idxTop(~act)-1;

      % Remap the dofs
      allDofs = 1:obj.ndofs;
      actDofs = obj.getActiveDof;
      [tf, loc] = ismember(allDofs, actDofs);
      dofs = loc;
      dofs(~tf) = 0;

      npoint = sum(4*idxTop+4,"all");
      data(1:npoint,1)=0.;

      count = 0;
      for col=1:prod(obj.ncells(1:2))
        idI = mod(col-1,obj.ncells(1))+1;
        idJ = (col-idI)/obj.ncells(1)+1;

        count = count+1;
        topLay = idxTop(col);
        acc = 0.;
        for lay = 1:topLay
          cell = dofs(obj.dof(idI,idJ,lay));
          acc = acc+comp(cell);
          data(4*count+1:4*count+4,1) = acc;
          count=count+1;
        end
      end
    end

    function [dof,pos] = getDofFromLay(obj,nlay)
      dof = obj.dof(:,:,nlay);
      map = ismember(dof,obj.inactDof);
      dof(map) = 0;
      map = dof(:)~=0;
      [dof,id] = sort(dof(map));
      pos = (1:length(dof))';
      pos = pos(map);
      pos = pos(id);
    end

    function dofs = getBord(obj,label)
      % GETBORDDOFS Returns dofs for the border.
      %
      % Outputs:
      %   id        - DOFs at the boundary
      %
      % Supported labels:
      %   'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'

      switch lower(label)
        case "xmin"
          dofs = obj.getBordXMin();
        case "xmax"
          dofs = obj.getBordXMax();
        case "ymin"
          dofs = obj.getBordYMin();
        case "ymax"
          dofs = obj.getBordYMax();
        case "zmin"
          dofs = obj.getBordZMin();
        case "zmax"
          dofs = obj.getBordZMax();
        otherwise
          dofs = [];
      end
    end

    function dofs = reScale(obj,dof)
      % Remap the dofs
      allDofs = 1:obj.ndofs;
      actDofs = obj.getActiveDof;
      [tf, loc] = ismember(allDofs, actDofs);
      ddofs = loc;
      ddofs(~tf) = 0;
      dofs = ddofs(dof);      
    end



    function dofs = getBordXMin(obj)
      mapH = reshape(obj.laysByCol,obj.ncells(1:2));
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
      mapH = reshape(obj.laysByCol,obj.ncells(1:2));
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
      mapH = reshape(obj.laysByCol,obj.ncells(1:2));
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
      mapH = reshape(obj.laysByCol,obj.ncells(1:2));
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
      % map = obj.laysByCol~=0;
      % cellID = sub2ind(obj.ncells,idI(map),idJ(map),idK(map));

      cellID = sub2ind(obj.ncells,idI,idJ,idK);
      dofs = obj.dof(cellID);
    end

    function dofs = getBordZMax(obj)
      % Find I and J positions.
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      % idK = findIdLaysVaryCell(obj);
      % idK = idK(:) - ~varCellAct;

      % layId = obj.laysByCol;
      % map = ismember(obj.dof,obj.inactDof);
      % pos = find(map);
      % [~,id] = sort(obj.dof(map));
      % [idI,idJ,idK]=ind2sub(obj.ncells,pos(id));

      % % Find K position.
      doftmp = obj.dof;
      doftmp(ismember(obj.dof,obj.inactDof))= 0;
      mask = doftmp ~= 0;
      idK = sum(mask, 3);
      idK = idK(:);

      % restrict to where exist a column.
      map = idK > 0;
      idI = idI(map);
      idJ = idJ(map);
      idK = idK(map);

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

    function out = getMapFromDofs(obj,dofs)
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
      minval = min(obj.laysByCol);
      if minval~=1
        loc = obj.laysByCol == minval;
        [idI, idJ] = obj.getIJLay();
        pos = sub2ind(obj.ncells, idI(loc), idJ(loc), obj.laysByCol(loc)-1);
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
      atTop = obj.laysByCol == obj.ncells(3);
      newlayer = any(and(atTop(:),map));
      if newlayer
        obj.ncells(3) = obj.ncells(3)+1;
        obj.dof(:,:,end+1) = zeros(obj.ncells(1:2));
      end

      % update the dof for the newest cells.
      dofs = obj.ndofs+1:obj.ndofs+cellsTadd;
      [idI,idJ] = obj.getIJLay;
      pos = sub2ind(obj.ncells,idI(map),idJ(map),obj.laysByCol(map)+1);
      obj.dof(pos) = dofs';
      obj.ndofs = obj.ndofs + cellsTadd;
      obj.laysByCol(map) = obj.laysByCol(map)+1;
      obj.inactDof(end+1:end+cellsTadd) = dofs;
    end

    function npts = getNumberPoints(obj)
      % GETNUMBERPOINTS Returns total number of grid points.
      npts = (obj.ncells(1)+1)*(obj.ncells(2)+1)*(obj.ncells(3)+1);
    end

  end


  methods (Access = private)

    function id = getBordX(obj,idx,nelm)
      id=zeros(nelm,1);

      count=1;
      for i=1:obj.ncells(1)
        for k=1:obj.ncells(3)
          pos = obj.ncells(1)*(i-1)+idx(i);
          if k>obj.laysByCol(pos)
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
          if k>obj.laysByCol(pos)
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

    function layId = findIdLaysVaryCell(obj)
      % Find the index of the layer in which the cell vary the height.
      % layId = obj.ncells(3)*ones(obj.ncells(1:2));
      layId = obj.laysByCol;
      map = ismember(obj.dof,obj.inactDof);
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [idI,idJ,idK]=ind2sub(obj.ncells,pos(id));
      cols = obj.ncells(1)*(idJ-1)+idI;
      for ln = 1:length(cols)
        val = idK(ln);
        if val < layId(cols(ln))
          layId(cols(ln)) = idK(col);
        end
      end
    end


  end
end