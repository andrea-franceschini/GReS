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

    function [id,faceArea] = getBordCell(obj,label)
      % GETBORDCELL Returns boundary cell information.
      %
      % Outputs:
      %   id        - Cell DOFs at the boundary
      %   faceArea  - Face areas
      %
      % Supported labels:
      %   'x0', 'xm', 'y0', 'ym', 'z0', 'zm'

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
          [id,faceArea] = getBordX(obj,J_idx,nelm);

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
          [id,faceArea] = getBordX(obj,J_idx,nelm);

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
          [id,faceArea] = getBordY(obj,I_idx,nelm);

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
          [id,faceArea] = getBordY(obj,I_idx,nelm);

        case "z0"
          [id,faceArea] = getBordZ(obj,true);

        case "zm"
          [id,faceArea] = getBordZ(obj,false);

        otherwise
          id = [];
          faceArea = [];
          return
      end
    end

    function ncells = getNumberCells(obj)
      % GETNUMBERCELLS Returns the number of active cells.
      ncells = sum(obj.dof~=0,"all");
      % ncells = sum(obj.columnsHeight, 'all');
    end

    function npts = getNumberPoints(obj)
      % GETNUMBERPOINTS Returns total number of grid points.
      npts = (obj.ncells(1)+1)*(obj.ncells(2)+1)*(obj.ncells(3)+1);
    end

    function coord = getCoordCenter(obj, cellIds)
      % GETCOORDCENTER Returns cell center coordinates.
      %
      % Input:
      %   cellIds - Active cell DOF indices
      %
      % Output:
      %   coord   - (x,y,z) coordinates of cell centers
      
      % Generate 3D matrices for each coordinate component
      [X, Y, Z] = ndgrid(obj.coordX(1:end-1)+diff(obj.coordX)/2., ...
        obj.coordY(1:end-1)+diff(obj.coordY)/2., ...
        obj.coordZ(1:end-1)+diff(obj.coordZ)/2.);

      % Extract the coordinates for the requested linear indices
      actCells = obj.dof ~= 0;
      coord = [X(actCells), Y(actCells), Z(actCells)];
      % dofs = obj.getActiveDofs();
      % coord = coord(dofs,:);
      % coord = coord(cellIds,:);
    end

    function vols = computeVols(obj,dofs,dl)
      % COMPUTEVOLS Computes volumes of all active cells.
      %
      % Notes:
      %   - Uses vectorized broadcasting
      %   - Assumes orthogonal grid

      if ~exist("dofs","var")
        dofs = obj.getActiveDofs;
      end
      if ~exist("dl","var")
        dl = zeros(length(dofs),1);
      end
      ind = find(ismember(obj.dof, dofs));
      
      dx = diff(obj.coordX(:));
      dy = reshape(diff(obj.coordY(:)), 1, []);
      dz = reshape(diff(obj.coordZ(:)), 1, 1, []);
      tmp = ones(length(dx),length(dy),length(dz));
      dz = dz.*tmp;

      vAll = dx.*dy.*tmp;
      vols = vAll(ind).*(dz(ind)+dl(dofs));
    end

    function ijk = getIJKfromCellID(obj,cellID)
      % GETIJKFROMCELLID Converts linear DOF to (i,j,k).
      ind = find(ismember(obj.dof, cellID));
      [i,j,k]=ind2sub(obj.ncells,ind);
      % [i,j,k]=ind2sub(obj.ncells,cellID);
      ijk=[i,j,k];
    end

    function cellID = getCellIDfromIJK(obj,i,j,k)
      % GETCELLIDFROMIJK Converts (i,j,k) to linear DOF.
      cellID = sub2ind(obj.ncells,i,j,k);
      % map = obj.dof(obj.dof~=0);
      % cellID=cellID(map);
    end

    function [idI,idJ,idK] = getIJKTop(obj)
      [idI, idJ] = obj.getIJLay;
      % idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      % idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = obj.columnsHeight;
    end

    function dofs = getActiveDofs(obj)
      % Returns a list of all non-zero Cell IDs.
      dofs = obj.dof(obj.dof ~=0);
    end

    function dofs = getDofsFromIJK(obj,ijk)
      if nargin == 1
        dofs = obj.dof(:);
      else
        pos = sub2ind(obj.ncells, ijk(:,1), ijk(:,2), ijk(:,3));
        dofs = obj.dof(pos);
      end
    end

    function neigh = getNeigh(obj,ijk)
      % GETNEIGH Returns neighboring cell DOFs.
      %
      % Order:
      %   [x-, x+, y-, y+, z-, z+]

      neigh = zeros(size(ijk,1),6);

      actTmp = ijk(:,1)~=1;
      pos = sub2ind(obj.ncells, ijk(actTmp,1)-1, ijk(actTmp,2), ijk(actTmp,3));
      neigh(actTmp,3) =  obj.dof(pos);
      actTmp = ijk(:,1)~=obj.ncells(1);
      pos = sub2ind(obj.ncells, ijk(actTmp,1)+1, ijk(actTmp,2), ijk(actTmp,3));
      neigh(actTmp,4) =  obj.dof(pos);

      actTmp = ijk(:,2)~=1;
      pos = sub2ind(obj.ncells, ijk(actTmp,1), ijk(actTmp,2)-1, ijk(actTmp,3));
      neigh(actTmp,1) =  obj.dof(pos);
      actTmp = ijk(:,2)~=obj.ncells(2);
      pos = sub2ind(obj.ncells, ijk(actTmp,1), ijk(actTmp,2)+1, ijk(actTmp,3));
      neigh(actTmp,2) =  obj.dof(pos);

      actTmp = ijk(:,3)~=1;
      pos = sub2ind(obj.ncells, ijk(actTmp,1), ijk(actTmp,2), ijk(actTmp,3)-1);
      neigh(actTmp,5) =  obj.dof(pos);
      actTmp = ijk(:,3)~=obj.ncells(3);
      pos = sub2ind(obj.ncells, ijk(actTmp,1), ijk(actTmp,2), ijk(actTmp,3)+1);
      neigh(actTmp,6) =  obj.dof(pos);
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
      pos = obj.getCellIDfromIJK(idI(map),idJ(map),obj.columnsHeight(map)+1);
      obj.dof(pos) = (obj.ndofs+1:obj.ndofs+cellsTadd)';
      obj.ndofs = obj.ndofs + cellsTadd;
      obj.columnsHeight(map) = obj.columnsHeight(map)+1;
    end

    function dims = getDims(obj,ijk)
      segmX = diff(obj.coordX);
      segmY = diff(obj.coordY);
      segmZ = diff(obj.coordZ);
      dims(:,1)=segmX(ijk(:,1));
      dims(:,2)=segmY(ijk(:,2));
      dims(:,3)=segmZ(ijk(:,3));
    end

    function conect = getConectByIJK(obj,idI,idJ,idK)
      % GETCONECTBYIJK Returns VTK hexahedral connectivity.

      ndivX = obj.ncells(1)+1;
      ndivXY = (obj.ncells(1)+1)*(obj.ncells(2)+1);

      conect=zeros(size(idI,1),8);

      % conect(:,1)=ndivXY*(idK-1)+4*(ndivX-1)*(idJ-1)+2*(idI-1);

      conect(:,1)=ndivXY*(idK-1)+ndivX*(idJ-1)+idI;
      conect(:,2)=ndivXY*(idK-1)+ndivX*(idJ-1)+idI+1;
      conect(:,3)=ndivXY*(idK-1)+ndivX*idJ+idI+1;
      conect(:,4)=ndivXY*(idK-1)+ndivX*idJ+idI;
      conect(:,5)=ndivXY*idK+ndivX*(idJ-1)+idI;
      conect(:,6)=ndivXY*idK+ndivX*(idJ-1)+idI+1;
      conect(:,7)=ndivXY*idK+ndivX*idJ+idI+1;
      conect(:,8)=ndivXY*idK+ndivX*idJ+idI;
    end

    function conect = getConectByIJKOk(obj,idI,idJ,idK)
      % GETCONECTBYIJK Returns VTK hexahedral connectivity.

      ndivX = obj.ncells(1)+1;
      ndivXY = (obj.ncells(1)+1)*(obj.ncells(2)+1);

      conect=zeros(size(idI,1),8);
      conect(:,1)=ndivXY*(idK-1)+ndivX*(idJ-1)+idI;
      conect(:,2)=ndivXY*(idK-1)+ndivX*(idJ-1)+idI+1;
      conect(:,3)=ndivXY*(idK-1)+ndivX*idJ+idI+1;
      conect(:,4)=ndivXY*(idK-1)+ndivX*idJ+idI;
      conect(:,5)=ndivXY*idK+ndivX*(idJ-1)+idI;
      conect(:,6)=ndivXY*idK+ndivX*(idJ-1)+idI+1;
      conect(:,7)=ndivXY*idK+ndivX*idJ+idI+1;
      conect(:,8)=ndivXY*idK+ndivX*idJ+idI;
    end

    function [coord, conect] = getMesh(obj)
      % GETMESH Returns mesh coordinates and connectivity.

      [XX, YY, ZZ] = ndgrid(obj.coordX, obj.coordY, obj.coordZ);
      coord = [XX(:), YY(:), ZZ(:)];

      dofs = obj.getActiveDofs;
      ijk = obj.getIJKfromCellID(dofs);
      conect = obj.getConectByIJK(ijk(:,1),ijk(:,2),ijk(:,3));
      % ddof = sub2ind(size(obj.dof), ijk(:,1), ijk(:,2), ijk(:,3));
      % ddof = obj.dof(ddof);
      % [~,idx]=sort(ddof);
      % conect=conect(idx,:);
    end

    function [idI, idJ, idK] = getIJK(obj)
      num=prod(obj.ncells);
      idI = reshape(repmat((1:obj.ncells(1))',prod(obj.ncells(2:3)),1),num,1);
      idJ = reshape(repmat(repelem((1:obj.ncells(2))',obj.ncells(1)),obj.ncells(3),1),num,1);
      idK = reshape(repelem((1:obj.ncells(3))',prod(obj.ncells(1:2))),num,1);
    end

    function dof = getMaxDofUnchanged(obj)
      minval = min(obj.columnsHeight);
      if minval~=1
        loc = obj.columnsHeight == minval;
        [idI, idJ] = obj.getIJLay();
        pos = sub2ind(obj.ncells, idI(loc), idJ(loc), obj.columnsHeight(loc)-1);
        dof = max(obj.dof(pos));
      else
        dof=0;
      end
    end

    function out = distMapOverDofsOld(obj,map)
      out = zeros(obj.ndofs,1);
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
      idK = obj.columnsHeight(:);
      dofTmp = obj.dof(sub2ind(obj.ncells,idI,idJ,idK));
      out(dofTmp) = map(:);
    end

    function out = distMapOverDofs(obj,map)
      out = zeros(obj.ndofs,1);
      for i=1:obj.ncells(1)
        for j=1:obj.ncells(2)
          dofTmp = obj.dof(sub2ind(obj.ncells,i,j,1:obj.ncells(3)));
          notZero =dofTmp~=0;
          out(dofTmp(notZero)) = map(i,j);
        end
      end
    end

    function dh = getCellHeight(obj)
      map = obj.dof ~= 0;
      pos = find(map);
      [~,id] = sort(obj.dof(map));
      [~,~,idK(id)]=ind2sub(obj.ncells,pos);
      segmZ = diff(obj.coordZ);
      dh=segmZ(idK);
    end

    function comp = getCompaction(obj,defByCell)
      comp = zeros(obj.ndofs,1);
      for i=1:obj.ncells(1)
        for j=1:obj.ncells(2)
          column = obj.dof(i,j,:);
          acc = 0.;
          for k=1:obj.ncells(3)
            dofId = column(k);
            if dofId ~= 0
              acc = acc + defByCell(dofId);
              comp(dofId)=acc;
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
        obj.gridClassic(data.Grid);
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

    function gridClassic(obj,data)
      % Internal initialization of grid, maps, and material layers.
      obj.ncells = str2num(data.division);
      
      % divisions = str2num(data.division);
      dim = str2num(data.size);
      obj.coordX = linspace(0,dim(1),obj.ncells(1)+1);
      obj.coordY = linspace(0,dim(2),obj.ncells(2)+1);
      obj.coordZ = linspace(0,dim(3),obj.ncells(3)+1);
    end

    function gridExplicit(obj,data)
      % Internal initialization of grid, maps, and material layers.
      if isfield(data,'xfile')
        if ~isfile(data.xfile)
          gridForSedimentation.errorGrid();
        end
        obj.coordX = load(data.xfile);
      else
        if ~isfield(data,'xcoord')
          gridForSedimentation.errorGrid();
        end
        obj.coordX = str2num(data.xcoord);
      end

      if isfield(data,'yfile')
        if ~isfile(data.yfile)
          gridForSedimentation.errorGrid();
        end
        obj.coordY = load(data.yfile);
      else
        if ~isfield(data,'ycoord')
          gridForSedimentation.errorGrid();
        end
        obj.coordY = str2num(data.ycoord);
      end

      if isfield(data,'zfile')
        if ~isfile(data.zfile)
          gridForSedimentation.errorGrid();
        end
        obj.coordZ = load(data.zfile);
      else
        if ~isfield(data,'zcoord')
          gridForSedimentation.errorGrid();
        end
        obj.coordZ = str2num(data.zcoord);
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

    function [id,area] = getBordX(obj,idx,nelm)
      areas = diff(obj.coordX).*diff(obj.coordZ)';
      id=zeros(nelm,1);
      area=zeros(nelm,1);

      count=1;
      for i=1:obj.ncells(1)
        for k=1:obj.ncells(3)
          pos = obj.ncells(1)*(i-1)+idx(i);
          if k>obj.columnsHeight(pos)
            break;
          end
          id(count)=obj.dof(i,idx(i),k);
          area(count)=areas(i,k);
          count=count+1;
        end
      end
    end

    function [id,area] = getBordY(obj,idx,nelm)
      areas = diff(obj.coordY).*diff(obj.coordZ)';
      id=zeros(nelm,1);
      area=zeros(nelm,1);

      count=1;
      for j=1:obj.ncells(2)
        for k=1:obj.ncells(3)
          pos = obj.ncells(1)*(idx(j)-1)+j;
          if k>obj.columnsHeight(pos)
            break;
          end
          id(count)=obj.dof(idx(j),j,k);
          area(count)=areas(idx(j),k);
          count=count+1;
        end
      end
    end

    function [id,area] = getBordZ(obj,bord)
      activeCells = obj.columnsHeight~=0;
      nelm = sum(activeCells);
      areas = reshape(diff(obj.coordY).*diff(obj.coordX)',nelm,1);

      [idI,idJ] = obj.getIJLay;
      if bord   % Case at the bord z0
        idK = ones(prod(obj.ncells(1:2)),1);
      else      % Case at the bord zm
        idK = obj.columnsHeight;
      end
      pos = getCellIDfromIJK(obj,idI(activeCells),idJ(activeCells), ...
          idK(activeCells));
      id=obj.dof(pos);
      area=areas(activeCells);
    end

    function [idI, idJ] = getIJLay(obj)
      idI = repmat((1:obj.ncells(1))',obj.ncells(2),1);
      idJ = repelem((1:obj.ncells(2))',obj.ncells(1));
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