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
            case 5 % 5 nodes
                a = 0.906179845938664;
                b = 0.5384693101056831;
                c = 0.0;
                wa = 0.2369268850561891;
                wb = 0.4786286704993665;
                wc = 0.568888888888888;
                obj.coord1D = [-a; -b; c; b; a];
                obj.weight1D = [wa; wb; wc; wb; wa];
            case 6 % 6 nodes
                a = 0.9324695142031521;
                b = 0.6612093864662645;
                c = 0.2386191860831969;
                wa = 0.1713244923791704;
                wb = 0.3607615730481386;
                wc = 0.4679139345726910;
                obj.coord1D = [-a; -b; -c; c; b; a];
                obj.weight1D = [wa; wb; wc; wc; wb; wa];
            case 15 % 15 nodes
                a = 0.9879925180204854;
                b = 0.9372733924007060;
                c = 0.8482065834104272;
                d = 0.7244177313601701;
                e = 0.5709721726085388;
                f = 0.3941513470775634;
                g = 0.2011940939974345;
                h = 0.0;
                wa = 0.030753241996117;
                wb = 0.0703660474881081;
                wc = 0.1071592204671719;
                wd = 0.1395706779261543;
                we = 0.1662692058169939;
                wf = 0.1861610000155622;
                wg = 0.1984314853271116;
                wh = 0.2025782419255613;
                obj.coord1D = [-a; -b; -c; -d; -e; -f; -g; h; g; f; e; d; c; b; a];
                obj.weight1D = [wa; wb; wc; wd; we; wf; wg; wh; wg; wf; we; wd; wc; wb; wa];

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