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
  end
  

  methods (Access = public)
    % Class constructor method
      function [outVar1,outVar2] = getDerBasisFAndDet(obj,el,flOut)   % mat,dJWeighed
      %       findJacAndDet(obj,el);  % OUTPUT: J and obj.detJ
      % Find the Jacobian matrix of the isoparametric map and its determinant
      %
      % Possible ways of calling this function are:
      %    1) [mat,dJWeighed] = getDerBasisFAndDet(obj,el,1)
      %    2) mat = getDerBasisFAndDet(obj,el,2)
      %    3) dJWeighed = getDerBasisFAndDet(obj,el,3)

      nodes = obj.grid.getCellNodes(el);
      coords = obj.grid.coordinates(nodes,:);
      
      if obj.GaussPts.nNode < 2
        % faster version for single GP rule
        mat = [ones(obj.nNode,1), coords];
        obj.detJ = det(mat);
        dJw = obj.GaussPts.weight*obj.detJ;
        invMat = inv(mat);
        N = invMat(2:obj.nNode,:);
      else
        [N, dJw] = mxGetDerBasisAndDet(obj.Jref,coords,obj.GaussPts.weight);
      end

      switch flOut
        case 1
          outVar1 = N;
          outVar2 = dJw';
        case 2
          outVar1 = N;
        case 3
          outVar1 = dJw';
      end
      if (flOut == 1 || flOut == 3) && obj.GaussPts.nNode > 1
        obj.detJ = (dJw./obj.GaussPts.weight)';
      end
      end

      function N = getBasisFinGPoints(obj)
        N = obj.Nref;
      end

      function computeProperties(obj)
        idTetra = find(obj.grid.cellVTKType == obj.vtkType);
        [vol,cellCent] = findVolumeAndCentroid(obj,idTetra);
        obj.grid.cellCentroid(idTetra,:) = cellCent;
        obj.grid.cellVolume(idTetra,:) = vol;
      end

      %   Elements volume calculation
      function [vol, cellCentroid] = findVolumeAndCentroid(obj, idTetra)

        % fix element orientation if some tetrahedra have wrong numbering

        if nargin == 1
          idTetra = find(obj.grid.cells.VTKtype == obj.vtkType);
        end

        tetraNodes = obj.getCellNodes(idTetra);  % [nTetra × 4]

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

      function volNod = findNodeVolume(obj,el)
        volNod = 0.25*obj.grid.cellVolume(el);
        volNod = repelem(volNod,obj.nNode);      
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

    function setElement(obj)
      obj.GaussPts = Gauss(obj.vtkType,obj.nGP);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
      findLocBasisF(obj);
      findLocDerBasisF(obj);
      FEM.setStrainMatrix(obj)
    end
  end
end