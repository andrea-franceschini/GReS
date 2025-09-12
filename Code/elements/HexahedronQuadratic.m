classdef HexahedronQuadratic < FEM
  % HEXAHEDRON element class
      %
      % NODE ORDERING ASSUMPTION (VTK, not available with msh grids)
      % Hexahedron (8 nodes):
      %
      % 4----11----3
      % |\         |\
      % |20    24  |19
      % 12 \  25   |10 \
      % |   8----15+---7
      % |21 |  27  | 22|
      % 1---+-9----2   |
      %  \  16    26\  14
      %  17 |  23    18|
      %    \|         \|
      %     5----13----6
        

    % map face
    % 1: 1-2-3-4
    % 2: 1-4-5-8
    % 3: 1-2-5-6
    % 4: 2-3-6-7
    % 5: 3-4-7-8
    % 6: 5-6-7-8

  properties (Constant)
    centroid = [0,0,0]

    coordLoc = [-1	-1	-1;
                 1	-1	-1;
                 1	 1	-1;
                -1   1	-1;
                -1	-1	 1;
                 1	-1	 1;
                 1	 1	 1;
                -1	 1	 1;
                 0	-1	-1;
                 1	 0	-1;
                 0	 1	-1;
                -1	 0	-1;
                 0	-1	 1;
                 1	 0	 1;
                 0	 1	 1;
                -1	 0	 1;
                -1	-1	 0;
                 1	-1	 0;
                 1	 1	 0;
                -1	 1	 0;
                -1	 0	 0;
                 1	 0	 0;
                 0	-1	 0;
                 0	 1	 0;
                 0	 0	-1;
                 0	 0	 1;
                 0	 0	 0];

    vtkType = 29
    nNode = 27
    nFace = 6
  end

  methods (Access = public)

    function [outVar1,outVar2] = getDerBasisFAndDet(obj,el,flOut)   % mat,dJWeighed
      %       findJacAndDet(obj,el);  % OUTPUT: J and obj.detJ
      % Find the Jacobian matrix of the isoparametric map and its determinant
      %
      % Possible ways of calling this function are:
      %    1) [mat,dJWeighed] = getDerBasisFAndDet(obj,el,1)
      %    2) mat = getDerBasisFAndDet(obj,el,2)
      %    3) dJWeighed = getDerBasisFAndDet(obj,el,3)

      coords = obj.mesh.coordinates(obj.mesh.cells(el,:),:);
      [N, dJw] = mxGetDerBasisAndDet(obj.Jref,coords,obj.GaussPts.weight);

      switch flOut
        case 1
          outVar1 = N;
          outVar2 = dJw';
        case 2
          outVar1 = N;
        case 3
          outVar1 = dJw';
      end
      if flOut == 1 || flOut == 3
        obj.detJ = (dJw./obj.GaussPts.weight)';
      end
    end

    function [outVar1,outVar2] = getDerBubbleBasisFAndDet(obj,el,flOut)   % mat,dJWeighed
      %       findJacAndDet(obj,el);  % OUTPUT: J and obj.detJ
      % Find the Jacobian matrix of the isoparametric map and its determinant
      %
      % Possible ways of calling this function are:
      %    1) [mat,dJWeighed] = getDerBasisFAndDet(obj,el,1)
      %    2) mat = getDerBasisFAndDet(obj,el,2)
      %    3) dJWeighed = getDerBasisFAndDet(obj,el,3)

      J = pagemtimes(obj.Jref,obj.mesh.coordinates(obj.mesh.cells(el,:),:));

      if flOut == 3 || flOut == 1
        %         obj.detJ = arrayfun(@(x) det(J(:,:,x)),1:obj.GaussPts.nNode);
        for i=1:obj.GaussPts.nNode
          obj.detJ(i) = det(J(:,:,i));
        end
        if flOut == 1
          outVar2 = obj.detJ.*(obj.GaussPts.weight)';
        elseif flOut == 3
          outVar1 = obj.detJ.*(obj.GaussPts.weight)';
        end
      end
      if flOut == 2 || flOut == 1
        for i=1:obj.GaussPts.nNode
          J(:,:,i) = inv(J(:,:,i));
        end
        outVar1 = pagemtimes(J,obj.Jb);
      end
    end


    function N1Mat = getBasisFinGPoints(obj) 
      N1Mat = obj.Nref;
    end

    function NbMat = getBubbleBasisFinGPoints(obj)
      NbMat = obj.Nb;
    end


    function [vol,cellCentroid] = findVolumeAndCentroid(obj,idHexa)
      % Find the volume of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      vol = zeros(length(idHexa),1);
      %       obj.volNod = zeros(obj.mesh.nNodes,1);
      cellCentroid = zeros(length(idHexa),3);
      i = 0;
      for el = idHexa'
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el,3);
        vol(i) = sum(dJWeighed);
        assert(vol(i)>0,'Volume less than 0');
        gPCoordinates = getGPointsLocation(obj,el);
        cellCentroid(i,:) = obj.detJ * gPCoordinates/vol(i);
      end
    end

    function nodeVol = findNodeVolume(obj,el)
      dJWeighed = obj.getDerBasisFAndDet(el,3);
      nodeVol = obj.Nref'*dJWeighed';
    end

    function gPCoordinates = getGPointsLocation(obj,el)
      % Get the location of the Gauss points in the element in the physical
      % space
      gPCoordinates = obj.Nref*obj.mesh.coordinates(obj.mesh.cells(el,:),:);
    end

    function computeProperties(obj)
      idHexa = find(obj.mesh.cellVTKType == obj.vtkType);
      [vol,cellCent] = findVolumeAndCentroid(obj,idHexa);
      obj.mesh.cellCentroid(idHexa,:) = cellCent;
      obj.mesh.cellVolume(idHexa,:) = vol;
    end

  end

  methods (Access = protected)
    function findLocDerBasisF(obj,varargin)
      % Compute derivatives in the reference space for all Gauss points
      if isempty(varargin)
        np = obj.GaussPts.nNode;
        c = obj.GaussPts.coord;
      else
        % compute basis at given reference point (xi,eta)
        c = varargin{1};
        np = size(c,1);
      end

      b = {@(x) 0.5*x*(x-1); @(x) 0.5*x*(x+1); @(x) 1-x^2};
      gb = {@(x) x-0.5; @(x) x+0.5; @(x) -2*x};
      idx = obj.tensProd(obj.coordLoc);

      obj.Jref = zeros(3,obj.nNode,np);

      for i = 1:obj.nNode

        [j,k,l] = deal(idx(i,1),idx(i,2),idx(i,3));
        obj.Jref(1,i,:) = arrayfun(@(g) ...
          gb{j}(c(g,1)).*b{k}(c(g,2)).*b{l}(c(g,3)),1:np);

        obj.Jref(2,i,:) = arrayfun(@(g) ...
          b{j}(c(g,1)).*gb{k}(c(g,2)).*b{l}(c(g,3)),1:np);

        obj.Jref(3,i,:) = arrayfun(@(g) ...
          b{j}(c(g,1)).*b{k}(c(g,2)).*gb{l}(c(g,3)),1:np);
      end
    end

    function findLocBasisF(obj,varargin)
      % Find the value the basis functions take at the Gauss points
      if isempty(varargin)
        np = obj.GaussPts.nNode;
        c = obj.GaussPts.coord;
      else
        % compute basis at given reference point (xi,eta)
        c = varargin{1};
        np = size(c,1);
      end
      
      b = {@(x) 0.5*x*(x-1); @(x) 0.5*x*(x+1); @(x) 1-x^2};
      N = @(i,j,k,x) b{i}(x(1)).*b{j}(x(2)).*b{k}(x(3));
      idx = obj.tensProd(obj.coordLoc);

      obj.Nref = zeros(np,obj.nNode);
      
      for i = 1:obj.nNode
        [j,k,l] = deal(idx(i,1),idx(i,2),idx(i,3)); 
        obj.Nref(:,i) = arrayfun(@(g) N(j,k,l,c(g,:)),1:np);
      end
    end

%     function findLocBubbleBasisF(obj)
%       % Find the value the basis functions take at the Gauss points
%       g = obj.GaussPts.coord;
%       bub = @(x,y) 1-g(x,y)^2;
%       val = @(x,y,z) 0.5 + 0.5*x*g(y,z);
% 
%       obj.Nb = zeros(obj.GaussPts.nNode,obj.nFace);
% 
%       obj.Nb(:,1) =  arrayfun(@(i) bub(i,1).*bub(i,2).*val(-1,i,3),...
%         1:obj.GaussPts.nNode);
%       obj.Nb(:,2) =  arrayfun(@(i) val(-1,i,1).*bub(i,2).*bub(i,3),...
%         1:obj.GaussPts.nNode);
%       obj.Nb(:,3) =  arrayfun(@(i) bub(i,1).*val(-1,i,2).*bub(i,3),...
%         1:obj.GaussPts.nNode);
%       obj.Nb(:,4) =  arrayfun(@(i) val(+1,i,1).*bub(i,2).*bub(i,3),...
%         1:obj.GaussPts.nNode);
%       obj.Nb(:,5) =  arrayfun(@(i) bub(i,1).*val(+1,i,2).*bub(i,3),...
%         1:obj.GaussPts.nNode);
%       obj.Nb(:,6) =  arrayfun(@(i) bub(i,1).*bub(i,2).*val(+1,i,3),...
%         1:obj.GaussPts.nNode);
%     end

%     function findLocDerBubbleBasisF(obj)
%       % Find the value the basis functions take at the Gauss points
%       zeros(obj.mesh.nDim,obj.mesh.cellNumVerts(1),obj.GaussPts.nNode);
%       g = obj.GaussPts.coord;
%       bub = @(x,y) 1.0-g(x,y)^2;
%       val = @(x,y,z) 0.5 + 0.5*x*g(y,z);
%       gradbub = @(x,y) -2.0*g(x,y);
% 
%       % prefer clarity over compactness
% 
%       obj.Jb = zeros(3,6,obj.GaussPts.nNode);
%       for i = 1:obj.GaussPts.nNode
%         % face 1
%         obj.Jb(1,1,i) =  gradbub(i,1).*bub(i,2).*val(-1,i,3);
%         obj.Jb(2,1,i) =  bub(i,1).*gradbub(i,2).*val(-1,i,3);
%         obj.Jb(3,1,i) =  bub(i,1).*bub(i,2).*(-0.5);
% 
%         % face 2
%         obj.Jb(1,2,i) =  (-0.5).*bub(i,2).*bub(i,3);
%         obj.Jb(2,2,i) =  val(-1,i,1).*gradbub(i,2).*bub(i,3);
%         obj.Jb(3,2,i) =  val(-1,i,1).*bub(i,2).*gradbub(i,3);
% 
%         % face 3
%         obj.Jb(1,3,i) =  gradbub(i,1).*val(-1,i,2).*bub(i,3);
%         obj.Jb(2,3,i) =  bub(i,1).*(-0.5).*bub(i,3);
%         obj.Jb(3,3,i) =  bub(i,1).*val(-1,i,2).*gradbub(i,3);
% 
%         % face 4
%         obj.Jb(1,4,i) =  0.5.*bub(i,2).*bub(i,3);
%         obj.Jb(2,4,i) =  val(+1,i,1).*gradbub(i,2).*bub(i,3);
%         obj.Jb(3,4,i) =  val(+1,i,1).*bub(i,2).*gradbub(i,3);
% 
%         % face 5
%         obj.Jb(1,5,i) =  gradbub(i,1).*val(+1,i,2).*bub(i,3);
%         obj.Jb(2,5,i) =  bub(i,1).*0.5.*bub(i,3);
%         obj.Jb(3,5,i) =  bub(i,1).*val(+1,i,2).*gradbub(i,3);
% 
%         % face 6
%         obj.Jb(1,6,i) =  gradbub(i,1).*bub(i,2).*val(+1,i,3);
%         obj.Jb(2,6,i) =  bub(i,1).*gradbub(i,2).*val(+1,i,3);
%         obj.Jb(3,6,i) =  bub(i,1).*bub(i,2).*0.5;
%       end
%     end

    function setElement(obj)
      obj.GaussPts = Gauss(obj.vtkType,obj.nGP);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
      findLocBasisF(obj);
      findLocDerBasisF(obj);
%       findLocBubbleBasisF(obj);
%       findLocDerBubbleBasisF(obj);
      FEM.setStrainMatrix(obj);
    end
  end

  methods (Static)
    function input = tensProd(input)
      % map local coordinate to tensor product indices
      input(input== 0) = 3;
      input(input == +1) = 2;
      input(input == -1) = 1;
    end
  end


end