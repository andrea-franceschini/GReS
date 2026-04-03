classdef (Abstract) FEM < handle

  properties 
    grid
    indB
    indBbubble
    GaussPts
    Jref      % basis function ref. gradient
    Nref      % basis function ref.
    Jb        % bubble basis function ref. gradient
    Nb        % bubble basis function ref.
    interpOrd
    nGP
    detJ
  end

  properties (Access = protected)
    topolMap      % id of cells with this shape
  end
  

  methods (Access = public)

    % Abstract class constructor
    function obj = FEM(ng,grid)
        obj.nGP = ng;
      if nargin > 1
        obj.grid = grid;
      end
      setElement(obj);
    end

    getBasisFinGPoints(obj)
    getDerBasisFAndDet(obj)
    computeProperties(obj)
  end

  methods (Access = protected)
    setElement(obj)
    findLocBasisF(obj)
    findLocDerBasisF(obj)

    function cellNodes = getCellNodes(id)
      % get a list of cell nodes from the grid topology
      % return a matrix of size nCells x nNode
      if isMixed(obj.grid)
        nId = obj.grid.cells.connectivity.getArray(id);
        cellNodes = (reshape(nId,4,[]))';
      else
        cellNodes = obj.grid.cells.connectivity.getArray(id,:);
      end
    end
  end

  methods (Static)
    function setStrainMatrix(elem)
      Nb = elem.nNode*elem.GaussPts.nNode;
      elem.indB = Poromechanics.setStrainMatIndex(Nb);
      %
      Nb = elem.nFace*elem.GaussPts.nNode;
      elem.indBbubble = Poromechanics.setStrainMatIndex(Nb);
    end


    function coords = getElementCoords(elem,id)
      nodes = elem.mesh.surfaces(id,1:elem.nNode);
      coords = elem.mesh.coordinates(nodes,:);
    end

  end
end
