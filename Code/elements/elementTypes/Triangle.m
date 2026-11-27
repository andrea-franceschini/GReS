classdef Triangle < FiniteElementType
  % TRIANGLE element class

  properties (Constant)
    centroid = [1/3,1/3]
    coordLoc = [0 0; 1 0; 0 1];
    vtkType = 5
    nNode = 3
    nFace = 1
    minGaussOrder = 1
  end



  methods (Access = public)

    function gradN = getDerBasisF(obj,el)
      % compute derivatives of the basis functions for element in real
      % space 
      % this work only in 2D!

     % inv([ones(3,1), in])
     coords = obj.getElementCoords(el);
     coords = coords(:,1:2);
     inv_A = inv([ones(3,1), coords]);
     gradN = inv_A(2:3,:);

    end

    function dN = computeDerBasisF(obj,varargin)
      dN = obj.Jref;
    end

    function [out,weights] = getDerBasisFAndDet(obj,in)
      if isscalar(in)
        in = obj.getElementCoords(in);
      end
      J = [in(2,:) - in(1,:); in(3,:) - in(1,:)];
      if size(J,2) == 2
        obj.detJ = abs(det(J));
      else
        obj.detJ = norm(cross(J(1,:),J(2,:)));
      end
      weights = obj.detJ.*obj.GaussPts.weight';
      if nargout == 1
        out = weights;
      else
        out = pinv(J)*[-1 1 0; -1 0 1];
      end
    end

    function N1Mat = getBasisFinGPoints(obj)
      N1Mat = obj.Nref;
    end

    function N = computeBasisF(obj,coordList)
      % Find the value the basis functions take at some  reference points defined in
      % a list
      if size(coordList,1) > 1
        N = zeros(size(coordList,1),3);
        N(:,1) = 1-coordList(:,1)-coordList(:,2);
        N(:,2) = coordList(:,1);
        N(:,3) = coordList(:,2);
      else
        N = zeros(1,3);
        N(1) = 1-coordList(1)-coordList(2);
        N(2) = coordList(1);
        N(3) = coordList(2);
      end
    end

    function Nb = computeBubbleBasisF(obj,coordList)
      Nb = arrayfun(@(i) (1-coordList(i,1)-coordList(i,2)).*coordList(i,1).*coordList(i,2),...
        (1:size(coordList,1)));
      Nb = Nb';
    end

    function gPCoordinates = getGPointsLocation(obj,in)
      % Get the location of the Gauss points in the element in the physical
      % space
   
      if isscalar(in) % input is element id
        coords = obj.getElementCoords(in);
        gPCoordinates = obj.Nref*coords;
      else % input is coordinates matrix
        gPCoordinates = obj.Nref*in;
      end
    end

    function [area,cellCentroid] = getSizeAndCentroid(obj,idTri)
      % Find the Area of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      area = zeros(length(idTri),1);
      cellCentroid = zeros(length(idTri),size(obj.grid.coordinates,2));
      i = 0;
      for el = idTri'
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el);
        area(i) = sum(dJWeighed);
        assert(area(i)>0,'Volume less than 0');
        coord = obj.getElementCoords(el);
        cellCentroid(i,:) = 1/3*(sum(coord,1));
      end
    end

    function n = computeNormal(obj,idTri,varargin)
      % compute normal vector of triangle
      assert(isscalar(idTri),'Input id must be a scalar positive integer')
      nodeCoord = obj.getElementCoords(idTri);
      if size(nodeCoord,2) == 2
        nodeCoord(:,3) = 0;
      end
      v1 = nodeCoord(1,:) - nodeCoord(2,:);
      v2 = nodeCoord(2,:) - nodeCoord(3,:);
      n = cross(v1',v2');
      n = n/norm(n);
    end

    function areaNod = getNodeInfluence(obj,in)

      if isscalar(in)
        area = obj.grid.cells.volume(in);
      else
        [~,weights] = obj.getDerBasisFAndDet(in);
        area = sum(weights);
      end
      areaNod = repelem(area/obj.nNode,obj.nNode,1);
    end

    % function n_a = computeAreaNod(obj,surfMsh)
    %   % compute area associated to each node of a surface mesh
    %   n_a = zeros(max(surfMsh.surfaces,[],'all'),1);
    %   for i = 1:length(surfMsh.surfaces)
    %     a = findAreaAndCentroid(obj,i);
    %     n_a(surfMsh.surfaces(i,:)) = n_a(surfMsh.surfaces(i,:)) + a/3;
    %   end
    %   n_a = n_a(unique(surfMsh.surfaces));
    % end

  end

  methods (Access = protected)

    function findLocBasisF(obj)
      % Find the value the basis functions take at the Gauss points
      Ntmp = zeros(obj.GaussPts.nNode,3);
      Ntmp(:,1) = 1-obj.GaussPts.coord(:,1)-obj.GaussPts.coord(:,2);
      Ntmp(:,2) = obj.GaussPts.coord(:,1);
      Ntmp(:,3) = obj.GaussPts.coord(:,2);

      obj.Nref = reshape(Ntmp,obj.GaussPts.nNode,3);

    end

    function findLocDerBasisF(obj, varargin)
      % Find the value the basis functions take at the Gauss points
      obj.Jref = [-1 1 0; -1 0 1];
    end

    function setElement(obj)
      findLocBasisF(obj);
      findLocDerBasisF(obj);
    end

  end

  methods (Static)

    function out = checkInRange(coord,tol)

      if nargin < 2
        tol = 1e-5;
      end
      xi  = coord(:,1);
      eta = coord(:,2);

      out = (xi >= -tol) & (eta >= -tol) & (xi + eta <= 1 + tol);
    end

    function area = computeArea(coords)
      assert(all(size(coords)==[3 2]), ['Wrong size for input matrix of vertex coordinates. ' ...
        ' Expected size is 3 x 2']);

      x1 = coords(1,1); y1 = coords(1,2);
      x2 = coords(2,1); y2 = coords(2,2);
      x3 = coords(3,1); y3 = coords(3,2);

      area = 0.5 * abs( x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2) );

    end
  end



end
