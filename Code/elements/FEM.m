classdef (Abstract) FEM < handle

  properties
    mesh
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

  methods (Access = public)

    % Abstract class constructor
    function obj = FEM(ng,mesh)
        obj.nGP = ng;
      if nargin > 1
        obj.mesh = mesh;
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
  end

  methods (Static)
    function setStrainMatrix(elem)
      Nb = elem.nNode*elem.GaussPts.nNode;
      elem.indB = Poromechanics.setStrainMatIndex(Nb);
      %
      Nb = elem.nFace*elem.GaussPts.nNode;
      elem.indBbubble = Poromechanics.setStrainMatIndex(Nb);
    end

    function out = checkInRange(elem,coord)
      % check if reference coordinate in input are inside the reference
      % element
      rangeXi = [min(elem.coordLoc(:,1)) max(elem.coordLoc(:,1))];
      rangeEta = [min(elem.coordLoc(:,2)) max(elem.coordLoc(:,2))];
      out = true;
      tol = 1e-4; % to avoid issues with uneven number of gp
      if (coord(1) < rangeXi(1)-tol || coord(1) > rangeXi(2)+tol)
        out = false;
      elseif (coord(2) < rangeEta(1)-tol || coord(2) > rangeEta(2)+tol)
        out = false;
      end
    end


    function coords = getElementCoords(elem,id)
      nodes = elem.mesh.surfaces(id,1:elem.nNode);
      coords = elem.mesh.coordinates(nodes,:);
    end

  end
end
