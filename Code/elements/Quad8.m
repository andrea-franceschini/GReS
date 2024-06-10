 classdef Quad8 < handle
  % QUADRILATERAL element class

  properties (Access = public)
%     vol
%     volNod
%     cellCentroid
  end
  
  properties (Access = private)   % PRIVATE
%
% NODE ORDERING ASSUMPTION (same as Gmsh output):
% Quadrilaterals (8 nodes):             
% 
%       | v
% 4-----7-----3            
% |     |     |           
% |     |     |           
% 8     +-----6---->                
% |           |   u  
% |           |      
% 1-----5-----2           
      

    coordLoc = [-1 -1;
                 1 -1;
                 1  1;
                -1  1;
                 0 -1;
                 1  0;
                 0  1;
                 -1 0]
    GaussPts
    J1
    mesh
%     vol
    J
    detJ
%     invJ;
    N1
  end

  methods (Access = public)
    % Class constructor method
    function obj = Quad8(msh,GPoints)
       % Calling the function to set element data
%        obj.setElementData(data);
       % ALLOCATE j
       obj.setQuad8(msh,GPoints);
%        findVolumeAndCentroid(obj);
%        obj.J = zeros(3,3,obj.GaussPts.nNode);
%        obj.detJ = zeros(obj.GaussPts.nNode,1);
    end
       
    function [outVar1,outVar2] = getDerBasisFAndDet(obj,el,flOut)   % mat,dJWeighed
%       findJacAndDet(obj,el);  % OUTPUT: J and detJ
      % Find the Jacobian matrix of the isoparametric map and the geometric
      % map
      % if nDim = 2, outVar2 is the norm of cross product (the determinant
      % would be 0)
      % if nDim = 3, outVar2 is the standard jacobian determinant
      % Possible ways of calling this function are:
      %    1) [mat,dJWeighed] = getDerBasisFAndDet(obj,el,1)
      %    2) mat = getDerBasisFAndDet(obj,el,2)
      %    3) dJWeighed = getDerBasisFAndDet(obj,el,3)
      obj.J = pagemtimes(obj.J1,obj.mesh.coordinates(obj.mesh.surfaces(el,:),1:obj.mesh.nDim));
      if flOut == 3 || flOut == 1
          %         obj.detJ = arrayfun(@(x) det(obj.J(:,:,x)),1:obj.GaussPts.nNode);
          if obj.mesh.nDim == 3
              for i=1:obj.GaussPts.nNode
                  obj.detJ(i) = norm(cross(obj.J(1,:,i),obj.J(2,:,i)),2);
              end
              elseif obj.mesh.nDim == 2
              for i=1:obj.GaussPts.nNode
                  % we still use detJ with abuse of notation 
                  obj.detJ(i) = det(obj.J(:,:,i));
              end
          end
          if flOut == 1
              outVar2 = obj.detJ.*(obj.GaussPts.weight)';
          elseif flOut == 3
              outVar1 = obj.detJ.*(obj.GaussPts.weight)';
          end
      end
      if flOut == 2 || flOut == 1
%         invJTmp = arrayfun(@(x) inv(obj.J(:,:,x)),1:obj.GaussPts.nNode,'UniformOutput',false);
%         obj.J = reshape(cell2mat(invJTmp),obj.mesh.nDim,obj.mesh.nDim,obj.GaussPts.nNode); %possibly we can overwrite J
%         clear invJTmp
        for i=1:obj.GaussPts.nNode
          obj.J(:,:,i) = inv(obj.J(:,:,i));
        end
        outVar1 = pagemtimes(obj.J,obj.J1);
      end
      end
   
    
%     function findJacAndDet(obj,el)
%       % Find the Jacobian matrix of the isoparametric map and its determinant
%       obj.J = pagemtimes(obj.J1,obj.mesh.coordinates(obj.mesh.cells(el,:),:));
%       obj.detJ = arrayfun(@(x) det(obj.J(:,:,x)),1:obj.GaussPts.nNode);
%     end
    
    
    function N1Mat = getBasisFinGPoints(obj)
      N1Mat = obj.N1;
    end

    
%     function mat = getDerBasisF(obj,el)
%       obj.J = pagemtimes(obj.J1,obj.mesh.coordinates(obj.mesh.cells(el,:),:));
%       invJTmp = arrayfun(@(x) inv(obj.J(:,:,x)),1:obj.GaussPts.nNode,'UniformOutput',false);
%       obj.invJ = reshape(cell2mat(invJTmp),obj.mesh.nDim,obj.mesh.nDim,obj.GaussPts.nNode); %possibly we can overwrite J
%       clear invJTmp
%       mat = pagemtimes(obj.invJ,obj.J1);
%     end
%     
%     function dJWeighed = getJacDet(obj)
%       obj.detJ = arrayfun(@(x) det(obj.J(:,:,x)),1:obj.GaussPts.nNode);
%       dJWeighed = obj.detJ.*(obj.GaussPts.weight)';
%     end
    
%     function v = getVolume(obj,el)
%       v = obj.vol(el);
%     end
    function [area,cellCentroid] = findAreaAndCentroid(obj,idHexa)
      % Find the Area of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      area = zeros(length(idHexa),1);
      cellCentroid = zeros(length(idHexa),3);
      i = 0;
      for el = idHexa
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el,3);
        area(i) = sum(dJWeighed);
        assert(area(i)>0,'Volume less than 0');
        gPCoordinates = getGPointsLocation(obj,el);
        cellCentroid(i,:) = obj.detJ * gPCoordinates/area(i);
      end
    end
    
    function nodeArea = findNodeArea(obj,idQuad)
      nodeArea = zeros(8*length(idQuad),1);
      ptr = 0;
      for el = idQuad
        dJWeighed = obj.getDerBasisFAndDet(el,3);
        nodeArea(ptr+1:ptr+8) = obj.N1'*dJWeighed';
        ptr = ptr + 8;
      end
    end

    function n = computeNormal(obj,idQuad)
        % compute normal vector of cell idQuad
        n = zeros(length(idQuad),3);
        for el = idQuad
            nodeCoord = obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
            v1 = nodeCoord(1,:) - nodeCoord(2,:);
            v2 = nodeCoord(2,:) - nodeCoord(3,:);
            n(el,:) = cross(v1,v2);
            n(el,:) = n(el,:)/norm(n(el,:));
        end
    end

    function gPCoordinates = getGPointsLocation(obj,el)
        % Get the location of the Gauss points in the element in the physical
        % space
        gPCoordinates = obj.N1*obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
    end

    function N = computeBasisF(obj, list)
        % Find the value the basis functions take at some  reference points defined in
        % a list
        Na = bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*list(i,1)).* ...
            (1+obj.coordLoc(j,2).*list(i,2)).* ...
            (obj.coordLoc(j,1).*list(i,1)...
            +obj.coordLoc(j,2).*list(i,2)-1), ...
            (1:size(list,1))',1:4);

        Nh = bsxfun(@(i,j) 1/2*(1-list(i,1).^2).*(1+obj.coordLoc(j,2).*list(i,2)), ...
            (1:size(list,1))',[5 7]);

        Nv = bsxfun(@(i,j) 1/2*(1-list(i,2).^2).*(1+obj.coordLoc(j,1).*list(i,1)),...
            (1:size(list,1))',[6 8]);

        if size(list,1)==1
            N = [Na; Nh(1); Nv(1); Nh(2);  Nv(2)]; 
        else
            N = [Na Nh(:,1) Nv(:,1) Nh(:,2) Nv(:,2)];
        end
    end

    function dN = computeDerBasisF(obj, list)
        % Compute derivatives in the reference space for all Gauss points
     % d(N)/d\csi
     d11 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
         (2*list(i,1).*obj.coordLoc(j,1)+...
         obj.coordLoc(j,2).*list(i,2)).*...
         (1+obj.coordLoc(j,2).*list(i,2)), ...
         (1:size(list,1))',1:4);

     d12 = bsxfun(@(i,j) 1/2*(-2*list(i,1)).*...
         (1+obj.coordLoc(j,2).*list(i,2)), ...
         (1:size(list,1))',[5 7]);

     d13 = bsxfun(@(i,j) 1/2*(1-list(i,2).^2).*...
         obj.coordLoc(j,1),(1:size(list,1))',[6 8]);

     %
     % d(N)/d\eta
     d21 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
         (2*list(i,2).*obj.coordLoc(j,2)+...
         obj.coordLoc(j,1).*list(i,1)).*...
         (1+obj.coordLoc(j,1).*list(i,1)), ...
         (1:size(list,1))',1:4);

     d22 = bsxfun(@(i,j) 1/2*(1-list(i,1).^2).*...
         obj.coordLoc(j,2),(1:size(list,1))',[5 7]);

     d23 = bsxfun(@(i,j) 1/2*(-2*list(i,2)).*...
         (1+obj.coordLoc(j,1).*list(i,1)),...
        (1:size(list,1))',[6 8]);

     %ordering columns
     if size(list,1)==1
         d1 = [d11; d12(1); d13(1); d12(2); d13(2)];
         d2 = [d21; d22(1); d23(1); d22(2); d23(2)];
     else
         d1 = [d11 d12(:,1) d13(:,1) d12(:,2) d13(:,2)];
         d2 = [d21 d22(:,1) d23(:,1) d22(:,2) d23(:,2)];
     end

     dN = [d1';d2'];
    end

  end

  methods (Access = private)
    function findLocDerBasisF(obj)
      % Compute derivatives in the reference space for all Gauss points
      obj.J1 = zeros(obj.mesh.nDim,obj.mesh.surfaceNumVerts(1),obj.GaussPts.nNode);
      %
      %d_1: vertices
      %d_2: horizontal edges
      %d_3: vertical edges

      % d(N)/d\csi
      d11 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
          (2*obj.GaussPts.coord(i,1).*obj.coordLoc(j,1)+...
          obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)).*...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
          (1:obj.GaussPts.nNode)',1:4);

      d12 = bsxfun(@(i,j) 1/2*(-2*obj.GaussPts.coord(i,1)).*...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
          (1:obj.GaussPts.nNode)',[5 7]);

      d13 = bsxfun(@(i,j) 1/2*(1-obj.GaussPts.coord(i,2).^2).*...
          obj.coordLoc(j,1),(1:obj.GaussPts.nNode)',[6 8]);

      %
      % d(N)/d\eta
      d21 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
          (2*obj.GaussPts.coord(i,2).*obj.coordLoc(j,2)+...
          obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).*...
          (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)), ...
          (1:obj.GaussPts.nNode)',1:4);

      d22 = bsxfun(@(i,j) 1/2*(1-obj.GaussPts.coord(i,1).^2).*...
          obj.coordLoc(j,2), (1:obj.GaussPts.nNode)',[5 7]);

      d23 = bsxfun(@(i,j) 1/2*(-2*obj.GaussPts.coord(i,2)).*...
          (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)),...
          (1:obj.GaussPts.nNode)',[6 8]);

      %ordering columns
      d1 = [d11 d12(:,1) d13(:,1) d12(:,2) d13(:,2)];
      d2 = [d21 d22(:,1) d23(:,1) d22(:,2) d23(:,2)];

      %
      obj.J1(1,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d1';
      obj.J1(2,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d2';
    end
    
    function findLocBasisF(obj, varargin)
      % Find the value the basis functions take at the Gauss points
      Na = bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).* ...
                     (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)).* ...
                     (obj.coordLoc(j,1).*obj.GaussPts.coord(i,1) + ...
                     obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)-1), ...
                     (1:obj.GaussPts.nNode)',1:4);
      Nh = bsxfun(@(i,j) 1/2*(1-obj.GaussPts.coord(i,1).^2).* ...
                     (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
                     (1:obj.GaussPts.nNode)',[5 7]);
      Nv = bsxfun(@(i,j) 1/2*(1-obj.GaussPts.coord(i,2).^2).* ...
                     (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)), ...
                     (1:obj.GaussPts.nNode)',[6 8]);
      % ordering columns
      obj.N1 = [Na Nh(:,1) Nv(:,1) Nh(:,2) Nv(:,2)];

      if obj.GaussPts.nNode == 1
          obj.N1 = obj.N1';
      end
    end
    
    function setQuad8(obj,msh,GPoints)
      obj.mesh = msh;
      obj.GaussPts = GPoints;
      findLocDerBasisF(obj);
      findLocBasisF(obj);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
    end
  end

  end
