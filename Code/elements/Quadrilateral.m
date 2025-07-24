classdef Quadrilateral < FEM
  % QUADRILATERAL element class
  %
  % NODE ORDERING ASSUMPTION (same as Gmsh output):
  % Quadrilatersl (4 nodes):
  %
  %       | v
  % 4-----|-----3
  % |     |     |
  % |     |     |
  % |     +-----|---->
  % |           |   u
  % |           |
  % 1-----------2

  properties (Constant)
    centroid = [0,0]
    coordLoc = [-1 -1;
      1 -1;
      1  1;
      -1  1;]
    vtkType = 9
    nNode = 4
    nFace = 1
  end


  methods (Access = public)

    function [outVar1,outVar2] = getDerBasisFAndDet(obj,in)   % mat,dJWeighed
      %       findJacAndDet(obj,el);  % OUTPUT: J and obj.detJ
      % way to call this method: if el is a scalar (element idx) the 3D
      % coordinates are retrieved by the corresponding mesh object. Only
      % the determinant is returned
      % if in is not scalar, it is a 4x2 list of 2D coordinates, the
      % gradient matrix and the determinant are returned

      if size(in,2) == 2
        % 2D - in is a given list of x-y coordinates for nodes
        J = pagemtimes(obj.Jref,in);
        for i=1:obj.GaussPts.nNode
          obj.detJ(i) = det(J(:,:,i));
          J(:,:,i) = inv(J(:,:,i));
        end
        outVar1 = pagemtimes(J,obj.Jref);
        outVar2 = obj.detJ.*(obj.GaussPts.weight)';
      else
        if isscalar(in)
          % 3D setting
          coord = FEM.getElementCoords(obj,in);
        elseif size(in,2)==3 % 3D list of coordinates in input
          coord = in;
        end
        J = pagemtimes(obj.Jref,coord);
        for i = 1:obj.GaussPts.nNode
          obj.detJ(i) = norm(cross(J(1,:,i),J(2,:,i)),2);
        end
        outVar1 = obj.detJ.*(obj.GaussPts.weight)';
      end
    end


    function N1Mat = getBasisFinGPoints(obj)
      N1Mat = obj.Nref;
    end

    function NbMat = getBubbleBasisFinGPoints(obj)
      NbMat = obj.Nb;
      NbMat = reshape(NbMat,obj.GaussPts.nNode,[]);
    end


    function [area,cellCentroid] = findAreaAndCentroid(obj,idQuad)
      % Find the Area of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      area = zeros(length(idQuad),1);
      cellCentroid = zeros(length(idQuad),3);
      i = 0;
      for el = idQuad'
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el);
        area(i) = sum(dJWeighed);
        assert(area(i)>0,'Volume less than 0');
        gPCoordinates = getGPointsLocation(obj,el);
        cellCentroid(i,:) = dJWeighed * gPCoordinates/area(i);
      end
    end



    function nodeArea = findNodeArea(obj,el)
        dJWeighed = obj.getDerBasisFAndDet(el);
        nodeArea = obj.Nref'*dJWeighed';
    end



    function n = computeNormal(obj,idQuad,pos)
      % compute normal vector of quadrilatral in specific reference point
      assert(isscalar(idQuad),'Input id must be a scalar positive integer')
      dN = obj.computeDerBasisF(pos);
      nodeCoord = obj.mesh.coordinates(obj.mesh.surfaces(idQuad,:),:);
      tang = dN*nodeCoord;
      crossTang = cross(tang(1,:)',tang(2,:)');
      n = crossTang/norm(crossTang); % normalize
    end



    function n_a = computeAreaNod(obj,surfMsh)
      % compute area associated to each node of a surface mesh
      n_a = zeros(max(surfMsh.surfaces,[],'all'),1);
      for i = 1:length(surfMsh.surfaces)
        n_a(surfMsh.surfaces(i,:)) = n_a(surfMsh.surfaces(i,:)) + findNodeArea(obj,i);
      end
      n_a = n_a(unique(surfMsh.surfaces));
    end


    function gPCoordinates = getGPointsLocation(obj,in)
      % Get the location of the Gauss points in the element in the physical
      % space
      if isscalar(in) % element id
        gPCoordinates = obj.Nref*FEM.getElementCoords(obj,in);
      else 
        assert(size(in,1)==4, ['List of coordinates in ' ...
          'input must be a 4x2 or 4x3 matrix'])
        gPCoordinates = obj.Nref*in;
      end
    end

    function N = computeBubbleBasisF(obj, coordList)
      % Find the value the bubble basis functions take at some  reference
      % points whose 2D coordinates are store in coord
      N = arrayfun(@(i) (1-coordList(i,1)^2).*(1-coordList(i,2)^2),(1:size(coordList,1)));
      N = N';
    end

    function computeProperties(obj)
      idQuad = find(obj.mesh.surfaceVTKType == obj.vtkType);
      [area,cellCent] = findAreaAndCentroid(obj,idQuad);
      obj.mesh.surfaceCentroid(idQuad,:) = cellCent;
      obj.mesh.surfaceArea(idQuad,:) = area;
    end
  end

  methods (Access = protected)
    function setElement(obj)
      obj.GaussPts = Gauss(obj.vtkType,obj.nGP);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
      findLocBasisF(obj);
      findLocDerBasisF(obj);
      findLocBubbleBasisF(obj);
    end

    function findLocDerBasisF(obj,varargin)
      % Compute derivatives in the reference space for all Gauss points
      obj.Jref = zeros(2,obj.nNode,obj.GaussPts.nNode);
      %
      % d(N)/d\csi
      if isempty(varargin)
        d1 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
          (1:obj.GaussPts.nNode)',1:obj.nNode);
        %
        % d(N)/d\eta
        d2 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
          (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)), ...
          (1:obj.GaussPts.nNode)',1:obj.nNode);
        % d2 = 1/8.*coord_loc(:,2).*(1+coord_loc(:,1).*pti_G(1)).*(1+coord_loc(:,3).*pti_G(3));
        %
        obj.Jref(1,1:obj.nNode,1:obj.GaussPts.nNode) = d1';
        obj.Jref(2,1:obj.nNode,1:obj.GaussPts.nNode) = d2';
      else
        refCoord = varargin{1};
        d1 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
          (1+obj.coordLoc(j,2).*refCoord(i,2)), ...
          (1:size(refCoord,1))',1:obj.nNode);
        %
        % d(N)/d\eta
        d2 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
          (1+obj.coordLoc(j,1).*refCoord(i,1)), ...
          (1:size(refCoord,1))',1:obj.nNode);
        % d2 = 1/8.*coord_loc(:,2).*(1+coord_loc(:,1).*pti_G(1)).*(1+coord_loc(:,3).*pti_G(3));
        %
        obj.Jref(1,1:obj.nNode,1:obj.GaussPts.nNode) = d1';
        obj.Jref(2,1:obj.nNode,1:obj.GaussPts.nNode) = d2';
      end
    end



    function findLocBasisF(obj, varargin)
      % Find the value the basis functions take at the Gauss points
      if isempty(varargin)
        obj.Nref = bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).* ...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
          (1:obj.GaussPts.nNode)',1:obj.nNode);
        if obj.GaussPts.nNode == 1
          obj.Nref = obj.Nref';
        end
      else
        % compute basis at given reference point (xi,eta)
        refCoord = varargin{1};
        bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*refCoord(i,1)).* ...
          (1+obj.coordLoc(j,2).*refCoord(i,2)), ...
          (1:size(refCoord,1))',1:nNOde);
      end
    end


    function findLocBubbleBasisF(obj)
      % Find the value the basis functions take at the Gauss points

      g = obj.GaussPts.coord;
      bub = @(x,y) 1-g(x,y)^2;

      obj.Nb = zeros(obj.GaussPts.nNode,1);
      obj.Nb =  arrayfun(@(i) bub(i,1).*bub(i,2),1:obj.GaussPts.nNode);
    end
  end

  methods (Static)
    function dN = computeDerBasisF(list)
      % Compute derivatives in the reference space for input list of
      % referencecoordinates
      % d(N)/d\csi

      dN = mxComputeDerBasisFQuad4(list);
%      cLoc = Quadrilateral.coordLoc;
%      nN = Quadrilateral.nNode;
% 
%      d1 = bsxfun(@(i,j) 1/4*cLoc(j,1).* ...
%        (1+cLoc(j,2).*list(i,2)), ...
%        (1:size(list,1)),1:nN);
%      %
%      % d(N)/d\eta
%      d2 = bsxfun(@(i,j) 1/4*cLoc(j,2).* ...
%        (1+cLoc(j,1).*list(i,1)), ...
%        (1:size(list,1)),1:nN);
%      %
%      dN = [d1';d2'];
    end

    function N = computeBasisF(coordList)
      % Find the value the basis functions take at some  reference points
      % whose 2D coordinates are store in coord
%       cLoc = Quadrilateral.coordLoc;
%       nN = Quadrilateral.nNode;

        % call to mex code
        N = mxComputeBasisFQuad4(coordList);

%       N = bsxfun(@(i,j) 1/4*(1+cLoc(j,1).*coordList(i,1)).* ...
%         (1+cLoc(j,2).*coordList(i,2)), ...
%         (1:size(coordList,1))',1:nN);
%       if size(N,2) ~= nN
%         N = N';
%       end
    end
  end

end
