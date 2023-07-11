classdef Gauss < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (Access = public)
    coord;
    weight;
    nNode;
  end
  
  properties (Access = private)
    nDim;
    nNode1D;
    cellType;   % According to VTK classification
    coord1D;
    weight1D;
  end
  
  methods (Access = public)
    function obj = Gauss(cType,n1D,nD)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      if nargin ~= 3
        error('Not enough input argument in Gauss class');
      end
      obj.setGauss(cType,n1D,nD);
      findGaussPoints(obj);
    end
  end
  
  methods (Access = private)
    
    function findGaussPoints(obj)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      switch obj.cellType
        case 10 % Tetrahedron
          error('Gauss points for tetrahedra not yet implemented.');
        case 12 % Hexahedra
          points1D(obj);
          switch obj.nDim
            case 1 % 1D
              obj.coord = obj.coord1D;
              obj.weight = obj.weight1D;
              obj.nNode = obj.nNode1D;
            case 2 % 2D
              [y, x] = meshgrid(obj.coord1D,obj.coord1D);
              obj.coord = [x(:), y(:)];
              obj.weight = bsxfun(@(a,b) a*b,obj.weight1D,(obj.weight1D)');
              obj.weight = obj.weight(:);
              obj.nNode = (obj.nNode1D)^2;
            case 3 % 3D
              [y, x, z] = meshgrid(obj.coord1D,obj.coord1D,obj.coord1D);
              obj.coord = [x(:), y(:), z(:)];
              weightGaussTmp = bsxfun(@(a,b) a*b,obj.weight1D,(obj.weight1D)');
              obj.weight = bsxfun(@(a,b) a*b,weightGaussTmp(:),(obj.weight1D)');
              obj.weight = obj.weight(:);
              obj.nNode = (obj.nNode1D)^3;
          end
        otherwise
          error('Gauss points for VTK element type %d not yet implemented.',obj.cellType);
      end
    end
    
    function points1D(obj)
      switch obj.nNode1D
        case 1 % 1 node
          obj.coord1D = 0.0;
          obj.weight1D = 2.0;
        case 2 % 2 nodes
          l = 1/sqrt(3);
          obj.coord1D = [-l; l];
          obj.weight1D = ones(2,1);
        case 3 % 3 nodes
          a = sqrt(0.6);
          b = 0.0;
          wa = 5/9;
          wb = 8/9;
          obj.coord1D = [-a; b; a];
          obj.weight1D = [wa; wb; wa];
        case 4 % 4 nodes
          a = 0.861136311594953;
          b = 0.339981043584856;
          wa = 0.347854845137454;
          wb = 0.652145154862546;
          obj.coord1D = [-a; -b; b; a];
          obj.weight1D = [wa; wb; wb; wa];
        otherwise
          error('Number of 1D nodes not supported yet');
      end
    end
    
    function setGauss(obj,cType,n1D,nD)
      obj.cellType = cType;
      obj.nNode1D = n1D;
      obj.nDim = nD;
    end
  end
end