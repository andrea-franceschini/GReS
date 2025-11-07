classdef EmbeddedFracture < handle
  % Class implementing Embedded finite element EFEM(0)
  % Reference: Cusini et al(2021).
  
  properties
    fractureArea
    mapFractureToCell
    fractureCenters
    n
    m1
    m2
    penaltyNormal
    penaltyTangential
    gN
    slip
    coes
    phi

  end
  
  methods
    function obj = EmbeddedFracture(inputStruct,mesh)
      obj.coes = inputStruct.Coulomb.coes;
      obj.phi = inputStruct.Coulomb.phi;
      processFractureElements(obj,inputStruct)

    end
%     
    function H = computeHeaviside(obj,fracId,nodeId)
      % 
    end

    function getDerTractionGap(obj)
      %
    end

    function processFractureElements(obj,inputStruct)
      plane = inputStruct.plane;
      obj.n = plane.norma;
      obj.m1 = plane.lengthVector;
      obj.m2 = plane.widthVector;
      dims = plane.size;

      % compute signed distance to spot cut cell
      % loop over cut cell and compute intersection of the plane with each
      % edge
      % check if at least one intersection point lies in the plane extent
      % if yes go on and order the point in ccw order sorting atan wrt
      % average center
      % with order list, gen area and centroid


  end
  end
end

