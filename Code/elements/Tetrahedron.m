classdef Tetrahedron < FEM
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
      
      if obj.GaussPts.nNode < 2
        % faster version 
        mat = [1 obj.mesh.coordinates(obj.mesh.cells(el,1),:);
                   1 obj.mesh.coordinates(obj.mesh.cells(el,2),:);
                   1 obj.mesh.coordinates(obj.mesh.cells(el,3),:);
                   1 obj.mesh.coordinates(obj.mesh.cells(el,4),:)];
        obj.detJ = det(mat);
        dJw = obj.GaussPts.weight*obj.detJ;
        invMat = inv(mat);
        N = invMat(2:obj.nNode,:);
      else
        coords = obj.mesh.coordinates(obj.mesh.cells(el,:),:);
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
        idTetra = find(obj.mesh.cellVTKType == obj.vtkType);
        [vol,cellCent] = findVolumeAndCentroid(obj,idTetra);
        obj.mesh.cellCentroid(idTetra,:) = cellCent;
        obj.mesh.cellVolume(idTetra,:) = vol;
      end

      %   Elements volume calculation
      function [vol,cellCentroid] = findVolumeAndCentroid(obj,idTetra)
        vol = zeros(length(idTetra),1);
        cellCentroid = zeros(length(idTetra),3);
        %       obj.volSign = ones(obj.mesh.nCells,1);
        %       obj.volNod = zeros(obj.mesh.nNodes,1);
        i = 0;
        for el = idTetra'
          i = i + 1;
          top = obj.mesh.cells(el,1:obj.nNode);
          coordMat = obj.mesh.coordinates(top,:);
          vol(i) = det([ones(obj.nNode,1) coordMat])/6;
          if vol(i) < 0
            %           obj.volSign(el) = -1;
            vol(i) = -vol(i);
          end
          cellCentroid(i,:) = sum(coordMat,1)/obj.nNode;
        end
      end

      function gPCoordinates = getGPointsLocation(obj,el)
        % Get the location of the Gauss points in the element in the physical
        % space
        gPCoordinates = obj.Nref*obj.mesh.coordinates(obj.mesh.cells(el,:),:);
      end

      function volNod = findNodeVolume(obj,el)
        volNod = 0.25*obj.mesh.cellVolume(el);
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
      FiniteElementLagrangian.setStrainMatrix(obj)
    end
  end
end