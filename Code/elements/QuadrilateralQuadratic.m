classdef QuadrilateralQuadratic < FEM
  % QUADRILATERAL element class
  %
  % NODE ORDERING ASSUMPTION (same as Gmsh output):
  % Quadrilatersl (4 nodes):
  %
  %       | v
  % 4-----7-----3
  % |  4  |  3  |
  % |     |     |
  % 8-----9-----6---->
  % |  1  |  2  |   u
  % |     |     |
  % 1-----5-----2

  properties (Constant)
    centroid = [0,0]
    coordLoc = [-1 -1;
                 1 -1;
                 1  1;
                -1  1;
                 0 -1;
                 1  0;
                 0  1;
                -1  0;
                 0  0];

    vtkType = 28
    nNode = 9
    nFace = 1

    nod2sub = [1,5,9,8;
               5,2,6,9;
               9,6,3,7;
               8,9,7,4];

    subMap =  [1  1;
              -1  1;
              -1 -1;
               1 -1;]
  end

  properties
    subQuad Quadrilateral
  end


  methods (Access = public)

    function [outVar1,outVar2] = getDerBasisFAndDet(obj,in)   % mat,dJWeighed
      %       findJacAndDet(obj,el);  % OUTPUT: J and obj.detJ
      % way to call this method: if el is a scalar (element idx) the 3D
      % coordinates are retrieved by the corresponding mesh object. Only
      % the determinant is returned
      % if in is not scalar, it is a 4x2 list of 2D coordinates, the
      % gradient matrix and the determinant are returned

      if isscalar(in)
        % 3D setting
        coord = FEM.getElementCoords(obj,in);
        J = pagemtimes(obj.Jref,coord);
        for i = 1:obj.GaussPts.nNode
          obj.detJ(i) = norm(cross(J(1,:,i),J(2,:,i)),2);
        end
        outVar1 = obj.detJ.*(obj.GaussPts.weight)';
      else
        % 2D - in is a given list of x-y coordinates for nodes
        J = pagemtimes(obj.Jref,in);
        for i=1:obj.GaussPts.nNode
          obj.detJ(i) = det(J(:,:,i));
          J(:,:,i) = inv(J(:,:,i));
        end
        outVar1 = pagemtimes(J,obj.Jref);
        outVar2 = obj.detJ.*(obj.GaussPts.weight)';
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



    function n = computeNormal(obj,idQuad,pos,idSub)
      % compute normal vector of quadrilatral in specific reference point
      assert(isscalar(idQuad),'Element input id must be a scalar positive integer')
      if nargin < 4
        dN = computeDerBasisF(obj,pos);
        nodeCoord = FEM.getElementCoords(obj,idQuad);
      else
        dN = Quadrilateral.computeDerBasisF(pos);
        nodeCoord = obj.getSubElementCoords(idQuad,idSub);
      end
      tang = dN*nodeCoord;
      crossTang = cross(tang(1,:)',tang(2,:)');
      n = crossTang/norm(crossTang);
    end



    function centroid = computeCentroid(obj,idQuad,idSub)
      assert(isscalar(idQuad),'Element input id must be a scalar positive integer')
      % compute centroid of specific subElement defined by idSub
      if nargin > 2
        nodeCoord = obj.getSubElementCoords(idQuad,idSub);
        dJWeighed = getDerBasisFAndDet(obj.subQuad,nodeCoord);
        gPCoordinates = getGPointsLocation(obj.subQuad,nodeCoord);
      elseif nargin == 2
        dJWeighed = getDerBasisFAndDet(obj,idQuad);
        gPCoordinates = getGPointsLocation(obj,idQuad);
      else
        error('Wrong number of input arguments.')
      end
      area = sum(dJWeighed);
      assert(area>0,'Volume less than 0');
      centroid = dJWeighed * gPCoordinates/area;
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
        assert(size(in,1)==4 && size(in,2)==2, ['List of coordinates in ' ...
          'input must be a 4x2 matrix'])
        gPCoordinates = obj.Nref*in;
      end
    end

    function subElement = createSubElement(obj)
      subElement = Quadrilateral(obj.nGP);
    end


    function N = computeBasisF(obj, coordList)
      % Find the value the basis functions take at some  reference points
      % whose 2D coordinates are store in coord
      % Find the value the basis functions take at the Gauss points
      N = mxComputeBasisFQuad9(coordList);
%       b1 = @(x) 0.5*x*(x-1);
%       b2 = @(x) 1-x^2;
%       b3 = @(x) 0.5*x*(x+1);
% 
%       c = coordList;
%       np = size(c,1);
%       N = zeros(np,obj.nNode);
%       N(:,1) = arrayfun(@(i) b1(c(i,1)).*b1(c(i,2)),1:np);
%       N(:,2) = arrayfun(@(i) b3(c(i,1)).*b1(c(i,2)),1:np);
%       N(:,3) = arrayfun(@(i) b3(c(i,1)).*b3(c(i,2)),1:np);
%       N(:,4) = arrayfun(@(i) b1(c(i,1)).*b3(c(i,2)),1:np);
%       N(:,5) = arrayfun(@(i) b2(c(i,1)).*b1(c(i,2)),1:np);
%       N(:,6) = arrayfun(@(i) b3(c(i,1)).*b2(c(i,2)),1:np);
%       N(:,7) = arrayfun(@(i) b2(c(i,1)).*b3(c(i,2)),1:np);
%       N(:,8) = arrayfun(@(i) b1(c(i,1)).*b2(c(i,2)),1:np);
%       N(:,9) = arrayfun(@(i) b2(c(i,1)).*b2(c(i,2)),1:np);
%       if np == 1
%         N = N';
%       end
    end

%     function N = computeBubbleBasisF(obj, coordList)
%       % Find the value the bubble basis functions take at some  reference
%       % points whose 2D coordinates are store in coord
%       N = arrayfun(@(i) (1-coordList(i,1)^2).*(1-coordList(i,2)^2),(1:size(coordList,1)));
%       N = N';
%     end

    function computeProperties(obj)
      idQuad = find(obj.mesh.surfaceVTKType == obj.vtkType);
      [area,cellCent] = findAreaAndCentroid(obj,idQuad);
      obj.mesh.surfaceCentroid(idQuad,:) = cellCent;
      obj.mesh.surfaceArea(idQuad,:) = area;
    end


    function dN = computeDerBasisF(obj, list)
      % Compute derivatives in the reference space for input list of
      % reference coordinates
      % d(N)/d\csi
      % 1D basis function
      dN = mxComputeDerBasisFQuad9(list);
%       b1 = @(x) 0.5*x*(x-1);
%       b2 = @(x) 1-x^2;
%       b3 = @(x) 0.5*x*(x+1);
%       % gradient of 1D basis functions
%       gb1 = @(x) 0.5*(2*x-1);
%       gb2 = @(x) -2*x;
%       gb3 = @(x) 0.5*(2*x+1);
% 
%       c = list;
%       np = size(c,1);
%       dN = zeros(2,obj.mesh.surfaceNumVerts(1),np);
%       dN(1,1,:) = arrayfun(@(i) gb1(c(i,1)).*b1(c(i,2)),1:np);
%       dN(2,1,:) = arrayfun(@(i) b1(c(i,1)).*gb1(c(i,2)),1:np);
% 
%       dN(1,2,:) = arrayfun(@(i) gb3(c(i,1)).*b1(c(i,2)),1:np);
%       dN(2,2,:) = arrayfun(@(i) b3(c(i,1)).*gb1(c(i,2)),1:np);
% 
%       dN(1,3,:) = arrayfun(@(i) gb3(c(i,1)).*b3(c(i,2)),1:np);
%       dN(2,3,:) = arrayfun(@(i) b3(c(i,1)).*gb3(c(i,2)),1:np);
% 
%       dN(1,4,:) = arrayfun(@(i) gb1(c(i,1)).*b3(c(i,2)),1:np);
%       dN(2,4,:) = arrayfun(@(i) b1(c(i,1)).*gb3(c(i,2)),1:np);
% 
%       dN(1,5,:) = arrayfun(@(i) gb2(c(i,1)).*b1(c(i,2)),1:np);
%       dN(2,5,:) = arrayfun(@(i) b2(c(i,1)).*gb1(c(i,2)),1:np);
% 
%       dN(1,6,:) = arrayfun(@(i) gb3(c(i,1)).*b2(c(i,2)),1:np);
%       dN(2,6,:) = arrayfun(@(i) b3(c(i,1)).*gb2(c(i,2)),1:np);
% 
%       dN(1,7,:) = arrayfun(@(i) gb2(c(i,1)).*b3(c(i,2)),1:np);
%       dN(2,7,:) = arrayfun(@(i) b2(c(i,1)).*gb3(c(i,2)),1:np);
% 
%       dN(1,8,:) = arrayfun(@(i) gb1(c(i,1)).*b2(c(i,2)),1:np);
%       dN(2,8,:) = arrayfun(@(i) b1(c(i,1)).*gb2(c(i,2)),1:np);
% 
%       dN(1,9,:) = arrayfun(@(i) gb2(c(i,1)).*b2(c(i,2)),1:np);
%       dN(2,9,:) = arrayfun(@(i) b2(c(i,1)).*gb2(c(i,2)),1:np);
    end

    function coord = getSubElementCoords(obj,idQuad,idSub)
      nodeSub = obj.mesh.surfaces(idQuad,obj.nod2sub(idSub,:));
      coord = obj.mesh.coordinates(nodeSub,:);
    end

    function xi_sub = mapref2sub(obj,xi_ref,sub)
      % map the reference coordinate of element to reference coordinate of
      % the subsegment space defined by 1D
      xi_sub = 2*xi_ref + obj.subMap(sub,:);
    end

    function xi_ref = mapsub2ref(obj,xi_sub,sub)
      % map reference coordinate of subsegment space defined by sub to
      % reference coordinate of element
      xi_ref = 0.5*(xi_sub - obj.subMap(sub,:));
    end
  end

  methods (Access = protected)
    function setElement(obj)
      obj.GaussPts = Gauss(obj.vtkType,obj.nGP);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
      findLocBasisF(obj);
      findLocDerBasisF(obj);
      obj.subQuad = Quadrilateral(obj.nGP);
%       findLocBubbleBasisF(obj);
    end

    function findLocDerBasisF(obj,varargin)
      
      % 1D basis function
      b1 = @(x) 0.5*x*(x-1);
      b2 = @(x) 1-x^2;
      b3 = @(x) 0.5*x*(x+1);
      % gradient of 1D basis functions
      gb1 = @(x) 0.5*(2*x-1);
      gb2 = @(x) -2*x;
      gb3 = @(x) 0.5*(2*x+1);

      if isempty(varargin)
        np = obj.GaussPts.nNode;
        c = obj.GaussPts.coord;
      else
        % compute basis at given reference point (xi,eta)
        c = varargin{1};
        np = size(c,1);
      end

      % Compute derivatives in the reference space for all Gauss points
      obj.Jref = zeros(2,obj.mesh.surfaceNumVerts(1),np);
      
      obj.Jref(1,1,:) = arrayfun(@(i) gb1(c(i,1)).*b1(c(i,2)),1:np);
      obj.Jref(2,1,:) = arrayfun(@(i) b1(c(i,1)).*gb1(c(i,2)),1:np);

      obj.Jref(1,2,:) = arrayfun(@(i) gb3(c(i,1)).*b1(c(i,2)),1:np);
      obj.Jref(2,2,:) = arrayfun(@(i) b3(c(i,1)).*gb1(c(i,2)),1:np);

      obj.Jref(1,3,:) = arrayfun(@(i) gb3(c(i,1)).*b3(c(i,2)),1:np);
      obj.Jref(2,3,:) = arrayfun(@(i) b3(c(i,1)).*gb3(c(i,2)),1:np);

      obj.Jref(1,4,:) = arrayfun(@(i) gb1(c(i,1)).*b3(c(i,2)),1:np);
      obj.Jref(2,4,:) = arrayfun(@(i) b1(c(i,1)).*gb3(c(i,2)),1:np);

      obj.Jref(1,5,:) = arrayfun(@(i) gb2(c(i,1)).*b1(c(i,2)),1:np);
      obj.Jref(2,5,:) = arrayfun(@(i) b2(c(i,1)).*gb1(c(i,2)),1:np);

      obj.Jref(1,6,:) = arrayfun(@(i) gb3(c(i,1)).*b2(c(i,2)),1:np);
      obj.Jref(2,6,:) = arrayfun(@(i) b3(c(i,1)).*gb2(c(i,2)),1:np);

      obj.Jref(1,7,:) = arrayfun(@(i) gb2(c(i,1)).*b3(c(i,2)),1:np);
      obj.Jref(2,7,:) = arrayfun(@(i) b2(c(i,1)).*gb3(c(i,2)),1:np);

      obj.Jref(1,8,:) = arrayfun(@(i) gb1(c(i,1)).*b2(c(i,2)),1:np);
      obj.Jref(2,8,:) = arrayfun(@(i) b1(c(i,1)).*gb2(c(i,2)),1:np);

      obj.Jref(1,9,:) = arrayfun(@(i) gb2(c(i,1)).*b2(c(i,2)),1:np);
      obj.Jref(2,9,:) = arrayfun(@(i) b2(c(i,1)).*gb2(c(i,2)),1:np);
    end



    function findLocBasisF(obj, varargin)
      % Find the value the basis functions take at the Gauss points
      b1 = @(x) 0.5*x*(x-1);
      b2 = @(x) 1-x^2;
      b3 = @(x) 0.5*x*(x+1);
      if isempty(varargin)
        np = obj.GaussPts.nNode;
        c = obj.GaussPts.coord;
      else
        % compute basis at given reference point (xi,eta)
        c = varargin{1};
        np = size(c,1);
      end
      obj.Nref = zeros(np,obj.nNode);
      obj.Nref(:,1) = arrayfun(@(i) b1(c(i,1)).*b1(c(i,2)),1:np);
      obj.Nref(:,2) = arrayfun(@(i) b3(c(i,1)).*b1(c(i,2)),1:np);
      obj.Nref(:,3) = arrayfun(@(i) b3(c(i,1)).*b3(c(i,2)),1:np);
      obj.Nref(:,4) = arrayfun(@(i) b1(c(i,1)).*b3(c(i,2)),1:np);
      obj.Nref(:,5) = arrayfun(@(i) b2(c(i,1)).*b1(c(i,2)),1:np);
      obj.Nref(:,6) = arrayfun(@(i) b3(c(i,1)).*b2(c(i,2)),1:np);
      obj.Nref(:,7) = arrayfun(@(i) b2(c(i,1)).*b3(c(i,2)),1:np);
      obj.Nref(:,8) = arrayfun(@(i) b1(c(i,1)).*b2(c(i,2)),1:np);
      obj.Nref(:,9) = arrayfun(@(i) b2(c(i,1)).*b2(c(i,2)),1:np);
      if np == 1
        obj.Nref = obj.Nref';
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

    function out = checkInRange(coord,tol)
      % check if reference coordinate in input are inside the reference
      % element

      if nargin < 2
        tol = 1e-5;
      end

      xi  = coord(:,1);
      eta = coord(:,2);

      out = (xi  >= -1 - tol) & (xi  <= 1 + tol) & ...
        (eta >= -1 - tol) & (eta <= 1 + tol);

    end

  end

end
