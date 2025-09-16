classdef Triangle < FEM
  % TRIANGLE element class

  properties (Constant)
    centroid = [1/3,1/3]
    coordLoc = [0 0; 1 0; 0 1];
    vtkType = 5
    nNode = 3
    nFace = 1
  end



  methods (Access = public)

    function [mat] = getDerBasisF(obj,el)
      % compute derivatives of the basis functions for element in real
      % space 
      % this work only in 2D!
      inv_A = inv([1 obj.mesh.coordinates(obj.mesh.surfaces(el,1),1:2);
        1 obj.mesh.coordinates(obj.mesh.surfaces(el,2),1:2);
        1 obj.mesh.coordinates(obj.mesh.surfaces(el,3),1:2)]);
      mat = inv_A(2:3,:);
    end

    function dN = computeDerBasisF(obj,varargin)
      dN = obj.Jref;
    end

     function varargout = getDerBasisFAndDet(obj,in)   % mat,dJWeighed
      % findJacAndDet(obj,el);  % OUTPUT: J and obj.detJ
      % way to call this method: if el is a scalar (element idx) the 3D
      % coordinates are retrieved by the corresponding mesh object. Only
      % the determinant is returned
      % if in is not scalar, it is a 4x2 list of 2 coordinates, the
      % gradient matrix and the determinant are returned

      if isscalar(in)
        % 3D setting
        assert(nargout==1,['Too many output argument. If input is scalar' ...
          ' 3D setting is assumed and only element jacobian is returned'])
        % jacobian is constant in a simplex
        coord = FEM.getElementCoords(obj,in);
        v1 = coord(1,:) - coord(2,:);
        v2 = coord(1,:) - coord(3,:);
        obj.detJ = norm(cross(v1,v2));
        varargout{1} = obj.detJ.*(obj.GaussPts.weight)';
      else
        % 2D setting: 'in' is a given list of x-y coordinates
        inv_A = inv([ones(3,1), in]);
        mat = inv_A(2:3,:);
        %v1 = norm(in(1,:)-in(2,:));
        %v2 = norm(in(1,:)-in(3,:));
        e1 = [in(1,:)-in(2,:), 0];
        e2 = [in(1,:)-in(3,:), 0];
        obj.detJ = norm(cross(e1',e2'));
        % jacobian is constant in a simplex
        if nargout == 2
          varargout{1} = mat;
          varargout{2} = obj.detJ.*(obj.GaussPts.weight)';
        else
          varargout{1} = obj.detJ*(obj.GaussPts.weight)';
        end
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
      if isscalar(in)
        gPCoordinates = obj.Nref*obj.mesh.coordinates(obj.mesh.surfaces(in,:),:);
      else
        gPCoordinates = obj.Nref*in;
      end
    end

    function [area,cellCentroid] = findAreaAndCentroid(obj,idTri)
      % Find the Area of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      area = zeros(length(idTri),1);
      cellCentroid = zeros(length(idTri),3);
      i = 0;
      for el = idTri'
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el);
        area(i) = sum(dJWeighed);
        assert(area(i)>0,'Volume less than 0');
        coord = FEM.getElementCoords(obj,el);
        cellCentroid(i,:) = 1/3*(sum(coord,1));
      end
    end

    function n = computeNormal(obj,idTri,varargin)
      % compute normal vector of triangle
      assert(isscalar(idTri),'Input id must be a scalar positive integer')
      nodeCoord = FEM.getElementCoords(obj,idTri);
      v1 = nodeCoord(1,:) - nodeCoord(2,:);
      v2 = nodeCoord(2,:) - nodeCoord(3,:);
      n = cross(v1',v2');
      n = n/norm(n);
    end

    function areaNod = findNodeArea(obj,el)
      areaNod = (1/obj.nNode)*obj.mesh.surfaceArea(el);
      areaNod = repelem(areaNod,obj.nNode,1);
    end

    function n_a = computeAreaNod(obj,surfMsh)
      % compute area associated to each node of a surface mesh
      n_a = zeros(max(surfMsh.surfaces,[],'all'),1);
      for i = 1:length(surfMsh.surfaces)
        a = findAreaAndCentroid(obj,i);
        n_a(surfMsh.surfaces(i,:)) = n_a(surfMsh.surfaces(i,:)) + a/3;
      end
      n_a = n_a(unique(surfMsh.surfaces));
    end


    function computeProperties(obj)
      % update the parent mesh object computing cell properties
      idTri = find(obj.mesh.surfaceVTKType == obj.vtkType);
      [area,centr] = findAreaAndCentroid(obj,idTri);
      obj.mesh.surfaceCentroid(idTri,:) = centr;
      obj.mesh.surfaceArea(idTri,:) = area;
    end
  end

  methods (Access = protected)

    function findLocBasisF(obj)
      % Find the value the basis functions take at the Gauss points
      Ntmp = zeros(obj.GaussPts.nNode,3);
      Ntmp(:,1) = 1-obj.GaussPts.coord(:,1)-obj.GaussPts.coord(:,2);
      Ntmp(:,2) = obj.GaussPts.coord(:,1);
      Ntmp(:,3) = obj.GaussPts.coord(:,2);
      obj.Nref = Ntmp;
      if obj.GaussPts.nNode == 1
        obj.Nref = obj.Nref';
      end
    end

    function findLocDerBasisF(obj, varargin)
      % Find the value the basis functions take at the Gauss points
      obj.Jref = [-1 1 0; -1 0 1];
    end

    function setElement(obj)
      obj.GaussPts = Gauss(obj.vtkType,obj.nGP);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
      findLocBasisF(obj);
      findLocDerBasisF(obj);
    end

  end
end
