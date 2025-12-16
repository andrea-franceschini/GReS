classdef structGrid
  %UNTITLED3 Summary of this class goes here
  %   Detailed explanation goes here

  properties
    X (:,1)
    Y (:,1)
    Z (:,1)
    centerX (:,1)
    centerY (:,1)
    centerZ (:,1)
  end

  properties (Access = private)
    npoints (1,3) int64
    ncells (1,3) int64
    ncellsSurf (1,3) int64 % XY, YZ, XZ
  end

  methods (Access = public)
    function obj = structGrid(varargin)
      %structGrid Construct an instance of this class
      %   Detailed explanation goes here
      if nargin == 0
        return
      end
      switch lower(varargin{1}.type)
        case "classic"
          obj=obj.constructorClassic(varargin{1});
        otherwise
      end
    end

    function [coord, conect] = getMesh(obj,varargin)
      if (nargin==1)
        coord = [];
        conect = [];
        return;
      end
      coord = obj.getMeshCoordColumn(varargin{1});
      conect = obj.getMeshConectColumn(varargin{1});
    end

    % Function to help this class
    function dim = getDimension(obj,varargin)
      % Return the size of the block
      if (nargin==1)
        dim = [obj.X(end),obj.Y(end),obj.Z(end)];
        return;
      end
      switch lower(varargin{1})
        case 'x'
          dim = obj.X(end);
        case 'y'
          dim = obj.Y(end);
        case 'z'
          dim = obj.Z(end);
        otherwise
          dim = [obj.X(end),obj.Y(end),obj.Z(end)];
      end
    end

    function npts = getNumberPoints(obj,varargin)
      % Function to return the size of the block
      if (nargin==1)
        npts = prod(obj.npoints);
        return;
      end
      switch lower(varargin{1})
        case 'x'
          npts = obj.npoints(1);
        case 'y'
          npts = obj.npoints(2);
        case 'z'
          npts = obj.npoints(3);
        otherwise
          npts = obj.npoints;
      end
    end

    function ncells = getNumberCells(obj,varargin)
      % Function to return the number of cells in the grid.
      if (nargin==1)
        ncells = obj.ncells;
        return;
      end
      switch lower(varargin{1})
        case 'x'
          ncells = obj.ncells(1);
        case 'y'
          ncells = obj.ncells(2);
        case 'z'
          ncells = obj.ncells(3);
        otherwise
          ncells = obj.ncells;
      end
    end

    function faces = getFacesAxis(obj,axis,layer)
      % start point in the face.
      pos = (layer-1)*obj.ncellsSurf(1)+1;
      pos = pos + (layer-1)*obj.npoints(1)*obj.ncells(2);
      pos = pos + (layer-1)*obj.npoints(2)*obj.ncells(1);
      switch axis
        case 'x'
          pos = pos + obj.ncellsSurf(1);
          nfacesRef = obj.npoints(1)*obj.ncells(2);
        case 'y'
          pos = pos + obj.ncellsSurf(1) + obj.npoints(1)*obj.ncells(2);
          nfacesRef = obj.npoints(2)*obj.ncells(1);
        case 'z'
          nfacesRef = obj.ncellsSurf(1);
      end
      faces = pos:(pos+nfacesRef-1);
    end

    function pts = getCellCenter(obj,axis)
      % % % % start point in the face.
      % % % pos = (layer-1)*obj.ncellsSurf(1)+1;
      % % % pos = pos + (layer-1)*obj.npoints(1)*obj.ncells(2);
      % % % pos = pos + (layer-1)*obj.npoints(2)*obj.ncells(1);
      switch axis
        case 'x'
          pts = obj.X(1:end-1)+diff(obj.X)/2;
        case 'y'
          pts = obj.Y(1:end-1)+diff(obj.Y)/2;
        case 'z'
          pts = obj.Z(1:end-1)+diff(obj.Z)/2;
      end
    end

    function surfs = getSurfs(obj,type,vecA,vecB)
      switch type
        case "x"
          segmA = diff(obj.Y);
          segmB = diff(obj.Z);
        case "y"
          segmA = diff(obj.X);
          segmB = diff(obj.Z);
        case "z"
          segmA = diff(obj.X);
          segmB = diff(obj.Y);
      end
      surfs = segmA(vecA).*segmB(vecB);
    end

    function vols = getVols(obj,posX,posY,posZ)
      segmX = diff(obj.X);
      segmY = diff(obj.Y);
      segmZ = diff(obj.Z);
      vols = segmX(posX).*segmY(posY).*segmZ(posZ);
    end

    % % % % function faces = getFacesLayer(obj)
    % % % %   % start point in the face.
    % % % %   pos = (layer-1)*obj.ncellsSurf(1)+1;
    % % % %   pos = pos + (layer-1)*obj.npoints(1)*obj.ncells(2);
    % % % %   pos = pos + (layer-1)*obj.npoints(2)*obj.ncells(1);
    % % % %   % 
    % % % %   % nelm = obj.ncellsSurf(1);
    % % % %   % z0face = pos:(pos+nelm-1);
    % % % %   % 
    % % % %   % pos = pos + obj.ncellsSurf(1);
    % % % %   % x0face = 
    % % % % 
    % % % % 
    % % % %   bot = obj.mesh.cells(:, [1, 4, 3, 2]);
    % % % %   top = obj.mesh.cells(:, [6, 7, 8, 5]);
    % % % %   est = obj.mesh.cells(:, [2, 3, 7, 6]);
    % % % %   wst = obj.mesh.cells(:, [5, 8, 4, 1]);
    % % % %   sth = obj.mesh.cells(:, [6, 5, 1, 2]);
    % % % %   nth = obj.mesh.cells(:, [3, 4, 8, 7]);
    % % % % end

    



  end

  methods (Access = private)
    function obj = constructorClassic(obj,data)
      obj.ncells = str2num(data.division);
      obj.npoints = obj.ncells+1;
      obj.ncellsSurf = [ obj.ncells(1)*obj.ncells(2), ...
        obj.ncells(2)*obj.ncells(3), obj.ncells(1)*obj.ncells(3)];
      % divisions = str2num(data.division);
      dim = str2num(data.size);
      obj.X = linspace(0,dim(1),obj.npoints(1));
      obj.Y = linspace(0,dim(2),obj.npoints(2));
      obj.Z = linspace(0,dim(3),obj.npoints(3));

      obj.centerX = obj.X(1:end-1)+diff(obj.X)/2;
      obj.centerY = obj.Y(1:end-1)+diff(obj.Y)/2;
      obj.centerZ = obj.Z(1:end-1)+diff(obj.Z)/2;
    end

    function coord = getMeshCoordColumn(obj,Kpos)
      % Return the coordinates from the layer kpos
      if Kpos>obj.npoints(3)
        coord = [];
        return
      end

      [XX, YY, ZZ] = ndgrid(obj.X, obj.Y, obj.Z(Kpos:end));
      coord = [XX(:), YY(:), ZZ(:)];
    end

    function conect = getMeshConectColumn(obj,Kpos)
      % Return the conectivities from the layer kpos
      % ndiv = [length(obj.arrayX),length(obj.arrayY),length(obj.arrayZ)]-1;
      % ndivX = length(obj.arrayX);
      % ndivXY = length(obj.arrayY)*ndivX;
      % ndivZ = length(obj.arrayZ)-1;
      if Kpos>obj.ncells(3)
        conect = [];
        return
      end

      ndivX = obj.npoints(1);
      ndivXY = prod(obj.npoints(1:2));

      % nelm = (ndiv(1)*ndiv(2))*(ndivZ-Kpos+1);
      nelm = prod(obj.ncells(1:2))*(obj.ncells(3)-Kpos+1);
      conect = zeros(nelm,8);
      pos=1;
      for k=Kpos:obj.ncells(3)
        for j=1:obj.ncells(2)
          for i=1:obj.ncells(1)
            conect(pos,1)=ndivXY*(k-1)+ndivX*(j-1)+i;
            conect(pos,2)=ndivXY*(k-1)+ndivX*(j-1)+i+1;
            conect(pos,3)=ndivXY*(k-1)+ndivX*j+i+1;
            conect(pos,4)=ndivXY*(k-1)+ndivX*j+i;
            conect(pos,5)=ndivXY*k+ndivX*(j-1)+i;
            conect(pos,6)=ndivXY*k+ndivX*(j-1)+i+1;
            conect(pos,7)=ndivXY*k+ndivX*j+i+1;
            conect(pos,8)=ndivXY*k+ndivX*j+i;
            pos=pos+1;
          end
        end
      end
    end

    

    function conect = getFaces2Cell(obj,Kpos)
      % Return the conectivities from the layer kpos
      if Kpos>ndivZ
        conect = [];
        return
      end


    end

    function data = getGeomData(obj,type,location)
      switch type
        case "volume"
        case ""
      end
    end

    

  end

end