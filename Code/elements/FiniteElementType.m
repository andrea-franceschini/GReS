classdef (Abstract) FiniteElementType < handle

  properties (Abstract,Constant)
    centroid
    coordLoc
    vtkType
    nNode
    nFace
    minGaussOrder
  end

  properties (Access=protected)
    grid
    indB              % indices of strain matrix
    indBbubble        % indices of bubble strain matrix
    GaussPts          % class for gauss point integration
    Jref              % basis function ref. gradient
    Nref              % basis function ref.
    Jb                % bubble basis function ref. gradient
    Nb                % bubble basis function ref.
    nGP               % gauss order
    detJ              % jacobian determinant
  end

  methods (Abstract,Access=public)

    getBasisFinGPoints(obj)       % compute basis functions
    getSizeAndCentroid(obj)
    getNodeInfluence(obj)         % compute area/volume influence to nodes

  end

  methods (Abstract,Access=protected)

    findLocBasisF(obj)
    findLocDerBasisF(obj)

  end


  methods (Access = public)

    % Abstract class constructor
    function obj = FiniteElementType(varargin)

      k = 1;

      if nargin > 0
        if isa(varargin{1},"Grid")
          obj.grid = varargin{1};
          k = 2;
        end
      end

      default = struct('gaussOrder',obj.minGaussOrder);
      g = readInput(default,varargin{k:end});

      if g.gaussOrder < 1
        g.gaussOrder = obj.minGaussOrder;
      end
      obj.GaussPts = Gauss(obj.vtkType,g.gaussOrder);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
      obj.setStrainMatrixIndex();

      setElement(obj);
    end


    function [N,dJw] = getDerBasisFAndDet(obj,coords)

      % coords: matrix of element coordinates of size nNode x 3

      % call mex file for basis function computation
      [N, dJw] = mxGetDerBasisAndDet(obj.Jref,coords,obj.GaussPts.weight);
      dJw = reshape(dJw,1,[]);

      obj.detJ = (dJw./obj.GaussPts.weight)';

    end

    
    function B = getStrainMatrix(obj,gradN)
      B = zeros(6,obj.nNode*obj.grid.nDim,obj.getNumbGaussPts);
      B(obj.indB(:,2)) = gradN(obj.indB(:,1));
    end


    function g = getGauss(obj)
      g = obj.GaussPts;
    end

    function nGP = getNumbGaussPts(obj)
      nGP = obj.GaussPts.nNode;
    end


    function coords = getElementCoords(obj,id)
      nodes = obj.grid.getSurfNodes(id);
      coords = obj.grid.coordinates(nodes,:);
    end

  end

  methods (Access = protected)

    function setElement(obj)
      findLocBasisF(obj);
      findLocDerBasisF(obj);
    end


    function setStrainMatrixIndex(obj)
      n = obj.nNode*obj.GaussPts.nNode;
      obj.indB = Poromechanics.setStrainMatIndex(n);
      % bubbles
      n = obj.nFace*obj.GaussPts.nNode;
      obj.indBbubble = Poromechanics.setStrainMatIndex(n);
    end


  end

  methods (Static)

    function elem = create(vtkId,varargin)
      % automatic dispatch of finite element tpyes
      switch vtkId
        case 5
          elem = Triangle(varargin{:});
        case 10
          elem = Tetrahedron(varargin{:});
        case 9
          elem = Quadrilateral(varargin{:});
        case 12
          elem = Hexahedron(varargin{:});
        case 28
          elem = QuadrilateralQuadratic(varargin{:});
        case 29
          elem = HexahedronQuadratic(varargin{:});
      end
    end




    function totGP = getTotGPinGrid(grid,ord,elemList)
      % return the total number of gauss points in the grid

      % gOrd: order of gauss integration (assumed constant for each element
      % type

      totGP = 0;

      for vtkId = grid.cells.vtkTypes

        elem = FiniteElementType.create(vtkId,grid,ord);
        ng = elem.getGauss.nNode;
        c = grid.getCellsByVTKId(vtkId);

        if nargin > 2
          c = intersect(c,elemList);
        end

        totGP = totGP + ng*numel(c);
      end
    end

  end

end
