classdef Tetrahedron < FiniteElementType
  % TETRAHEDRON element class

    properties (Constant)
    centroid = [0.25,0.25,0.25]
    coordLoc = [0 0 0;
      1 0 0;
      0 1 0;
      0 0 1];
    vtkType = 10
    nNode = 4
    nFace = 3
    minGaussOrder = 1
  end
  

  methods (Access = public)
    % Class constructor method


      function N = getBasisFinGPoints(obj)
        N = obj.Nref;
      end


      %   Elements volume calculation
      function [vol, cellCentroid] = getSizeAndCentroid(obj, idTetra)

        % fix element orientation if some tetrahedra have wrong numbering

        if nargin == 1
          idTetra = find(obj.grid.cells.VTKType == obj.vtkType);
        end

        tetraNodes = obj.grid.getCellNodes(idTetra);  % [nTetra × 4]

        X = obj.grid.coordinates(tetraNodes(:,1),:); % [nTetra × 3]
        Y = obj.grid.coordinates(tetraNodes(:,2),:);
        Z = obj.grid.coordinates(tetraNodes(:,3),:);
        W = obj.grid.coordinates(tetraNodes(:,4),:);

        cellCentroid = (X + Y + Z + W)/obj.nNode;

        v1 = Y - X;
        v2 = Z - X;
        v3 = W - X;

        vol = dot(v1, cross(v2,v3,2),2) / 6;

        degen = abs(vol) < 1e-10;

        if any(degen)
          error('Found %i degenerate tetrahedrons',sum(degen));
        end

        neg = vol < 0;

        % flip nodes 3 and 4 and change volume sign
        tetraNodes(neg,[3 4]) = tetraNodes(neg,[4 3]);

        if any(neg)
          gresLog().warning(1,"Found %i tetrahedra with negative determinant. Node ordering has been automatically fixed",sum(neg))
        end

        obj.grid.setConnectivity("cells",idTetra,tetraNodes);

      end

      function gPCoordinates = getGPointsLocation(obj,el)
        % Get the location of the Gauss points in the element in the physical
        % space
        nodes = obj.grid.getCellNodes(el);
        coords = obj.grid.coordinates(nodes,:);
        gPCoordinates = obj.Nref*coords;
      end

      function volNod = getNodeInfluence(obj,in)

        volNod = repelem(obj.grid.cells.volume(in)/obj.nNode,obj.nNode,1);
      end
      %     end
  end

  methods (Access = protected)
    function findLocBasisF(obj)
      g = obj.GaussPts.coord;
      ng = obj.GaussPts.nNode;
      obj.Nref = zeros(ng,obj.nNode);
      obj.Nref(:,1) = arrayfun(@(i) 1-g(i,1)-g(i,2)-g(i,3),1:ng);
      obj.Nref(:,2) = arrayfun(@(i) g(i,1),1:ng);
      obj.Nref(:,3) = arrayfun(@(i) g(i,2),1:ng);
      obj.Nref(:,4) = arrayfun(@(i) g(i,3),1:ng);
    end

    function findLocDerBasisF(obj)
      obj.Jref = [-1 1 0 0;
        -1 0 1 0;
        -1 0 0 1];
      obj.Jref = repmat(obj.Jref,1,1,obj.GaussPts.nNode);
    end
    % 
    % function setElement(obj)
    %   obj.GaussPts = Gauss(obj.vtkType,obj.nGP);
    %   obj.detJ = zeros(1,obj.GaussPts.nNode);
    %   findLocBasisF(obj);
    %   findLocDerBasisF(obj);
    %   FEM.setStrainMatrix(obj)
    % end
  end
end