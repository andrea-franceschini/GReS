classdef gridForSedimentation < handle
  % GRIDFORSEDIMENTATION
  % ------------------------------------------------------------------
  % Manages a structured 3D grid with dynamic vertical growth for
  % sedimentation-driven simulations.
  %
  % Responsibilities:
  %   - Store grid geometry and topology
  %   - Track active cells and column heights
  %   - Support dynamic addition of cells in the Z direction
  %
  % Notes:
  %   - Columns may have different heights
  %   - DOF numbering changes as the grid grows
  % ------------------------------------------------------------------

  properties
    % Grid geometry
    coordX (:,1)                 % X coordinates of grid nodes
    coordY (:,1)                 % Y coordinates of grid nodes
    coordZ (:,1)                 % Z coordinates of grid nodes
    ncells (1,3)                 % Number of cells in (x,y,z)

    % Dynamic grid control
    columnsHeight (:,:)          % Active height per column
    dof (:,:,:)                  % Cell DOF map (i,j,k) -> DOF ID
    ndofs                        % Number of active DOFs
  end

  methods (Access = public)
    function obj = gridForSedimentation(varargin)
      % GRIDFORSEDIMENTATION Constructor.
      %
      % Initializes the grid using XML input and material labels.
      %
      % Supported inputs:
      %   "XML"          - XML grid definition

      flagEscape=true;
      for k = 1:2:nargin
        switch lower(varargin{k})
          case "xml"
            data = varargin{k+1};
            flagEscape = false;
        end
      end

      flgInput = gridForSedimentation.checkInput(data);
      if flgInput
        gresLog().error("Grid for the sedimentation model not well defined!");
      end

      % Check if the class is gonna be constructed.
      if or(flagEscape,flgInput), return; end

      obj.constructor(data);
    end

    % DofManager related functions
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

        case "xmax"
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

        case "ymin"
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

        case "ymax"
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

        case "zmin"
          idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
          idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
          idK = ones(prod(obj.ncells(1:2)),1);

          % For case where some column do not exist
          % map = obj.columnsHeight~=0;
          % cellID = sub2ind(obj.ncells,idI(map),idJ(map),idK(map));

          cellID = sub2ind(obj.ncells,idI,idJ,idK);
          dofs = obj.dof(cellID);
        case "zmax"
          dofs = getTopDofs(obj);

        otherwise
          dofs = [];
          return
      end
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

    function dofs = getTopDofs(obj)  %<--- ok
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = obj.columnsHeight;

      % For case where some column do not exist
      % map = obj.columnsHeight~=0;
      % cellID = sub2ind(obj.ncells,idI(map),idJ(map),idK(map));

      cellID = sub2ind(obj.ncells,idI,idJ,idK);
      dofs = obj.dof(cellID);
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

    % Coordanates related functions
    function [dx,dy,dz] = getCellsDims(obj,dofs)
      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      segmX = diff(obj.coordX);
      segmY = diff(obj.coordY);
      segmZ = diff(obj.coordZ);
      map = ismember(obj.dof,dofs);
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [idI(id),idJ(id),idK(id)]=ind2sub(obj.ncells,pos);
      dx=segmX(idI');
      dy=segmY(idJ');
      dz=segmZ(idK');
    end

    function [x,y,z] = getCoordCenter(obj, dofs)
      % GETCOORDCENTER Returns cell center coordinates.
      %
      % Input:
      %   cellIds - Active cell DOF indices
      %
      % Output:
      %   coord   - (x,y,z) coordinates of cell centers

      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      
      % Generate 3D matrices for each coordinate component
      [x, y, z] = ndgrid(obj.coordX(1:end-1)+diff(obj.coordX)/2., ...
        obj.coordY(1:end-1)+diff(obj.coordY)/2., ...
        obj.coordZ(1:end-1)+diff(obj.coordZ)/2.);

      % Extract the coordinates for the requested linear indices
      map = ismember(obj.dof,dofs);
      id = sort(obj.dof(map));
      x = x(id);
      y = y(id);
      z = z(id);
    end

    function z = getCoordBottom(obj, dofs)
      % GETCOORDCENTER Returns cell center coordinates.
      %
      % Input:
      %   cellIds - Active cell DOF indices
      %
      % Output:
      %   coord   - (x,y,z) coordinates of cell centers

      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      
      % Generate 3D matrices for each coordinate component
      [~, ~, z] = ndgrid(obj.coordX(2:end),obj.coordY(2:end),obj.coordZ(2:end));

      % Generate 3D matrices for each coordinate component
      % [~, ~, z] = ndgrid(obj.coordX(2:end)+diff(obj.coordX)/2., ...
      %   obj.coordY(2:end)+diff(obj.coordY)/2., ...
      %   obj.coordZ(2:end)+diff(obj.coordZ)/2.);

      % Extract the coordinates for the requested linear indices
      map = ismember(obj.dof,dofs);
      id = sort(obj.dof(map));
      z = z(id);
    end

    function z = getCoordTop(obj, dofs)
      % GETCOORDCENTER Returns cell center coordinates.
      %
      % Input:
      %   cellIds - Active cell DOF indices
      %
      % Output:
      %   coord   - (x,y,z) coordinates of cell centers

      if ~exist("dofs","var")
        dofs = 1:obj.ndofs;
      end
      
      % Generate 3D matrices for each coordinate component
      [x, y, z] = ndgrid(obj.coordX(1:end-1), ...
        obj.coordY(1:end-1), obj.coordZ(1:end-1));

      % Extract the coordinates for the requested linear indices
      map = ismember(obj.dof,dofs);
      id = sort(obj.dof(map));
      z = z(id);
    end

    function zcol = getColumnMaxHeight(obj,dofs)
      if ~exist("dofs","var")
        dofs = (1:obj.ndofs)';
      end
      ndofs_eval = length(dofs);
      zcol = zeros(ndofs_eval,1);
      ind = find(ismember(obj.dof, dofs));
      [idI,idJ,~]=ind2sub(obj.ncells,ind);
      colMaxHei = reshape(obj.columnsHeight,obj.ncells(1:2));
      for i=1:ndofs_eval
        pos = colMaxHei(idI(i),idJ(i))+1;
        zcol(i) = obj.coordZ(pos,1);
      end
    end

    % Mesh related functions
    function [coord, conect] = getMesh(obj,dofs)
      if ~exist("dofs","var")
        dofs = (1:obj.ndofs)';
      end

      map = ismember(obj.dof,dofs);
      ref = find(map);
      [~,idx]=sort(obj.dof(map));
      [idI,idJ,idK]=ind2sub(obj.ncells,ref(idx));
      refdof = sub2ind(obj.ncells,idI,idJ,idK);
      dofs = obj.dof(refdof);

      conect = zeros(length(dofs),8);
      for i=1:8
        conect(:,i)=8*(dofs-1)+i;
      end

      coord = zeros(8*length(dofs),3);
      ddof = (1:length(dofs))';
      dofId = 8*(ddof-1);
      coord(dofId+1,1) = obj.coordX(idI(ddof)+0);
      coord(dofId+1,2) = obj.coordY(idJ(ddof)+0);
      coord(dofId+1,3) = obj.coordZ(idK(ddof)+0);

      coord(dofId+2,1) = obj.coordX(idI(ddof)+1);
      coord(dofId+2,2) = obj.coordY(idJ(ddof)+0);
      coord(dofId+2,3) = obj.coordZ(idK(ddof)+0);

      coord(dofId+3,1) = obj.coordX(idI(ddof)+1);
      coord(dofId+3,2) = obj.coordY(idJ(ddof)+1);
      coord(dofId+3,3) = obj.coordZ(idK(ddof)+0);

      coord(dofId+4,1) = obj.coordX(idI(ddof)+0);
      coord(dofId+4,2) = obj.coordY(idJ(ddof)+1);
      coord(dofId+4,3) = obj.coordZ(idK(ddof)+0);

      coord(dofId+5,1) = obj.coordX(idI(ddof)+0);
      coord(dofId+5,2) = obj.coordY(idJ(ddof)+0);
      coord(dofId+5,3) = obj.coordZ(idK(ddof)+1);

      coord(dofId+6,1) = obj.coordX(idI(ddof)+1);
      coord(dofId+6,2) = obj.coordY(idJ(ddof)+0);
      coord(dofId+6,3) = obj.coordZ(idK(ddof)+1);

      coord(dofId+7,1) = obj.coordX(idI(ddof)+1);
      coord(dofId+7,2) = obj.coordY(idJ(ddof)+1);
      coord(dofId+7,3) = obj.coordZ(idK(ddof)+1);

      coord(dofId+8,1) = obj.coordX(idI(ddof)+0);
      coord(dofId+8,2) = obj.coordY(idJ(ddof)+1);
      coord(dofId+8,3) = obj.coordZ(idK(ddof)+1);
    end

    function newlayer = grow(obj,map,height)
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
        obj.coordZ(end+1) = obj.coordZ(end)+height;
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

    % Other types
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
    function constructor(obj,data)
      % CONSTRUCTOR Internal grid initialization.
      %
      % Initializes:
      %   - Coordinates
      %   - DOF mapping
      %   - Column heights

      % Constructing the grid
      if strcmp(data.Grid.type,"classic")
        obj.ncells = data.Grid.division;
        dim = data.Grid.size;
        obj.coordX = linspace(0,dim(1),obj.ncells(1)+1);
        obj.coordY = linspace(0,dim(2),obj.ncells(2)+1);
        obj.coordZ = linspace(0,dim(3),obj.ncells(3)+1);
      elseif strcmp(data.Grid.type,"explicit")
        obj.gridExplicit(data.Grid);
      else
        gridForSedimentation.errorGrid();
      end

      % Initialaze some variables
      obj.dof = zeros(obj.ncells);
      obj.columnsHeight = zeros(prod(obj.ncells(1:2)),1);

      % Constructing the initial grid
      obj.initConf(data.Initial);

      obj.ndofs = sum(obj.dof~=0,"all");
    end

    function gridExplicit(obj,data)
      % Internal initialization of grid, maps, and material layers.
      if isfield(data,'xfile')
        obj.coordX = load(data.xfile);
      else
        obj.coordX = data.xcoord;
      end

      if isfield(data,'yfile')
        obj.coordY = load(data.yfile);
      else
        obj.coordY = data.ycoord;
      end

      if isfield(data,'zfile')
        obj.coordZ = load(data.zfile);
      else
        obj.coordZ = data.zcoord;
      end

      obj.ncells = [length(obj.coordX),length(obj.coordY),length(obj.coordZ)]-1;
    end

    function initConf(obj,data)
      % Defining the dof for the initial configuration
      switch data.type
        case "all"
          numcells = prod(obj.ncells);          
          obj.dof(1:numcells) = 1:numcells;
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
    end

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

    % Create the ijk position in the from 
    function [idI, idJ] = getIJLay(obj)
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
    end

    function [idI, idJ, idK] = getIJKfromDofs(obj,dofs)
      map = ismember(obj.dof,dofs);
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [idI,idJ,idK]=ind2sub(obj.ncells,pos(id));
    end

  end

  methods (Static)
    function flag = checkInput(input)
      % CHECKINPUT Validates grid and initial condition input.

      flag = false;
      if ~(or(isfield(input,'Grid'),isfield(input,'grid')))
        flag = true;
        gresLog().error("The grid parameters for the simulation is not defined!");
      end

      if ~or(isfield(input,'Initial'),isfield(input,'initial'))
          flag = true;
          gresLog().error("The material for each cell at the begin of the simulation is not defined!");
      end

    end

    function errorGrid()
      gresLog().error("The grid parameters for the simulation is not well defined!");
    end

  end

end