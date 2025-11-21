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


    function coords = getElementCoords(elem,id)
      nodes = elem.mesh.surfaces(id,1:elem.nNode);
      coords = elem.mesh.coordinates(nodes,:);
    end

  end
end
