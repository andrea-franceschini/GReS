classdef Triangle < handle
   % TRIANGLE element class

   properties (Access = private)
      mesh
      coordLoc = [-1 -1;
         1 -1;
         1  1;
         -1  1;]
      GaussPts
      N1
      J1
      detJ
   end

   properties (Access = public)
      cellCentroid;
   end

   methods (Access = public)
      % Class constructor method
      function obj = Triangle(mesh,varargin)
         % Calling the function to set element data
         if isempty(varargin{1})
            obj.setTriangle(mesh);
         else
            obj.setTriangle(mesh,varargin{1});
         end
         %       computeCellCentroid(obj);
      end
      %
      function [mat] = getDerBasisF(obj,el)
         % compute derivatives of the basis functions for element in real
         % space
         inv_A = inv([1 obj.mesh.coordinates(obj.mesh.surfaces(el,1),1:2);
            1 obj.mesh.coordinates(obj.mesh.surfaces(el,2),1:2);
            1 obj.mesh.coordinates(obj.mesh.surfaces(el,3),1:2)]);
         mat = inv_A(2:3,:);
      end
 
      function dN = computeDerBasisF(obj,varargin)
         dN = obj.J1;
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
         c = obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
         J = [c(2,1)-c(1,1) c(3,1)-c(1,1);
            c(2,2)-c(1,2) c(3,2)-c(1,2);
            c(2,3)-c(1,3) c(3,3)-c(1,3);];
         obj.detJ = norm(cross(J(:,1),J(:,2)),2);
         inv_A = inv([1 obj.mesh.coordinates(obj.mesh.surfaces(el,1),1:2);
            1 obj.mesh.coordinates(obj.mesh.surfaces(el,2),1:2);
            1 obj.mesh.coordinates(obj.mesh.surfaces(el,3),1:2)]);
         mat = inv_A(2:3,:);
         switch flOut
            case 1
               outVar1 = mat;
               outVar2 = obj.detJ*obj.GaussPts.weight';
            case 2
               outVar1 = mat;
            case 3
               outVar1 = obj.detJ*obj.GaussPts.weight';
         end
      end

      function N1Mat = getBasisFinGPoints(obj)
         N1Mat = obj.N1;
      end

      function N = computeBasisF(obj, list)
         % Find the value the basis functions take at some  reference points defined in
         % a list
         if size(list,1) > 1
            N = zeros(size(list,1),3);
            N(:,1) = 1-list(:,1)-list(:,2);
            N(:,2) = list(:,1);
            N(:,3) = list(:,2);
         else
            N = zeros(1,3);
            N(1) = 1-list(1)-list(2);
            N(2) = list(1);
            N(3) = list(2);
         end
      end

      function gPCoordinates = getGPointsLocation(obj,el)
         % Get the location of the Gauss points in the element in the physical
         % space
         gPCoordinates = obj.N1*obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
      end

      %   Elements volume calculation
      function vol = findVolume(obj,idTetra)
         vol = zeros(length(idTetra),1);
         %       obj.volSign = ones(obj.mesh.nCells,1);
         %       obj.volNod = zeros(obj.mesh.nNodes,1);
         i = 0;
         for el = idTetra
            i = i + 1;
            top = obj.mesh.surfaces(el,1:3);
            vol(i) = 0.5*det([1 obj.mesh.coordinates(top(1),1:2);
               1 obj.mesh.coordinates(top(2),1:2);
               1 obj.mesh.coordinates(top(3),1:2);]);
            if vol(i) < 0
               %           obj.volSign(el) = -1;
               vol(i) = -vol(i);
            end
            %         for i=1:obj.mesh.cellNumVerts(el)
            %           obj.volNod(top(i)) = obj.volNod(top(i)) + obj.vol(el)/obj.mesh.cellNumVerts(el);
            %         end
         end
         % Although it has no for loop, the following solution takes more
         % time!
         %       obj.vol = arrayfun(@(e) det([1 obj.mesh.coordinates(obj.mesh.cells(e,1),:);
         %                                    1 obj.mesh.coordinates(obj.mesh.cells(e,2),:);
         %                                    1 obj.mesh.coordinates(obj.mesh.cells(e,3),:);
         %                                    1 obj.mesh.coordinates(obj.mesh.cells(e,4),:)])/6,1:obj.mesh.nCells);
      end


      function [area,cellCentroid] = findAreaAndCentroid(obj,idTri)
         % Find the Area of the cells using the determinant of the Jacobian
         % of the isoparameric transformation
         area = zeros(length(idTri),1);
         cellCentroid = zeros(length(idTri),3);
         i = 0;
         for el = idTri
            i = i + 1;
            dJWeighed = getDerBasisFAndDet(obj,el,3);
            area(i) = sum(dJWeighed);
            assert(area(i)>0,'Volume less than 0');
            coord = obj.mesh.coordinates(obj.mesh.surfaces(idTri,:),:);
            cellCentroid(i,:) = 1/3*(sum(coord));
         end
      end

      function n = computeNormal(obj,idTri)
         % compute normal vector of cell idQuad
         n = zeros(length(idTri),3);
         for el = idTri
            nodeCoord = obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
            v1 = nodeCoord(1,:) - nodeCoord(2,:);
            v2 = nodeCoord(2,:) - nodeCoord(3,:);
            n(el,:) = cross(v1,v2);
            n(el,:) = n(el,:)/norm(n(el,:));
         end
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
   end

   methods (Access = private)
      function setTriangle(obj,mesh,varargin)
         obj.mesh = mesh;
         if ~isempty(varargin)
            obj.GaussPts = varargin{1};
            findLocBasisF(obj);
            findLocDerBasisF(obj);
         end
      end

      function findLocBasisF(obj, varargin)
         % Find the value the basis functions take at the Gauss points
         N1mat = zeros(obj.GaussPts.nNode,3);
         N1mat(:,1) = 1-obj.GaussPts.coord(:,1)-obj.GaussPts.coord(:,2);
         N1mat(:,2) = obj.GaussPts.coord(:,1);
         N1mat(:,3) = obj.GaussPts.coord(:,2);
         obj.N1 = N1mat;
         if obj.GaussPts.nNode == 1
            obj.N1 = obj.N1';
         end
      end

      function findLocDerBasisF(obj, varargin)
         % Find the value the basis functions take at the Gauss points
         obj.J1 = [-1 1 0; -1 0 1];
      end

      function computeCellCentroid(obj)
         obj.cellCentroid = sparse(repelem(1:obj.mesh.nSurfaces,obj.mesh.surfNumVerts), ...
            nonzeros((obj.mesh.surfaces)'),repelem((obj.mesh.surfNumVerts).^(-1), ...
            obj.mesh.surfNumVerts),obj.mesh.nSurfaces,obj.mesh.nNodes) ...
            * obj.mesh.coordinates;
      end
   end
end