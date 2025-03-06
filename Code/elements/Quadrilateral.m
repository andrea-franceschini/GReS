classdef Quadrilateral < handle
   % QUADRILATERAL element class


   properties (Access = private)   % PRIVATE
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

      coordLoc = [-1 -1;
         1 -1;
         1  1;
         -1  1;]
      GaussPts
      J1
      mesh
      J
      detJ
      N1
   end

   methods (Access = public)
      % Class constructor method
      function obj = Quadrilateral(msh,GPoints)
         % Calling the function to set element data
         %        obj.setElementData(data);
         % ALLOCATE j
         obj.setQuad(msh,GPoints);
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
         nodeArea = zeros(4*length(idQuad),1);
         ptr = 0;
         for el = idQuad
            dJWeighed = obj.getDerBasisFAndDet(el,3);
            nodeArea(ptr+1:ptr+4) = obj.N1'*dJWeighed';
            ptr = ptr + 4;
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

      function n_a = computeAreaNod(obj,surfMsh)
         % compute area associated to each node of a surface mesh
         n_a = zeros(max(surfMsh.surfaces,[],'all'),1);
         for i = 1:length(surfMsh.surfaces)
            n_a(surfMsh.surfaces(i,:)) = n_a(surfMsh.surfaces(i,:)) + findNodeArea(obj,i);
         end
         n_a = n_a(unique(surfMsh.surfaces));
      end

      function gPCoordinates = getGPointsLocation(obj,el)
         % Get the location of the Gauss points in the element in the physical
         % space
         gPCoordinates = obj.N1*obj.mesh.coordinates(obj.mesh.surfaces(el,:),:);
      end

      function N = computeBasisF(obj, list)
         % Find the value the basis functions take at some  reference points defined in
         % a list
         N = bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*list(i,1)).* ...
            (1+obj.coordLoc(j,2).*list(i,2)), ...
            (1:size(list,1))',1:obj.mesh.surfaceNumVerts(1));
         if size(N,2) ~= obj.mesh.surfaceNumVerts(1)
            N = N';
         end
      end

      function dN = computeDerBasisF(obj, list)
         % Compute derivatives in the reference space for all Gauss points
         % d(N)/d\csi
         d1 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
            (1+obj.coordLoc(j,2).*list(i,2)), ...
            (1:size(list,1)),1:obj.mesh.surfaceNumVerts(1));
         %
         % d(N)/d\eta
         d2 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
            (1+obj.coordLoc(j,1).*list(i,1)), ...
            (1:size(list,1)),1:obj.mesh.surfaceNumVerts(1));
         %
         dN = [d1';d2'];
      end

   end

   methods (Access = private)
      function findLocDerBasisF(obj)
         % Compute derivatives in the reference space for all Gauss points
         obj.J1 = zeros(obj.mesh.nDim,obj.mesh.surfaceNumVerts(1),obj.GaussPts.nNode);
         %
         % d(N)/d\csi
         d1 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,1).* ...
            (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
            (1:obj.GaussPts.nNode)',1:obj.mesh.surfaceNumVerts(1));
         %
         % d(N)/d\eta
         d2 = bsxfun(@(i,j) 1/4*obj.coordLoc(j,2).* ...
            (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)), ...
            (1:obj.GaussPts.nNode)',1:obj.mesh.surfaceNumVerts(1));
         % d2 = 1/8.*coord_loc(:,2).*(1+coord_loc(:,1).*pti_G(1)).*(1+coord_loc(:,3).*pti_G(3));
         %
         obj.J1(1,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d1';
         obj.J1(2,1:obj.mesh.surfaceNumVerts(1),1:obj.GaussPts.nNode) = d2';
         % the third row will remain 0 (shape functions are constant w.r.t the
         % third coordinate
      end

      function findLocBasisF(obj, varargin)
         % Find the value the basis functions take at the Gauss points
         obj.N1 = bsxfun(@(i,j) 1/4*(1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).* ...
            (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
            (1:obj.GaussPts.nNode)',1:obj.mesh.surfaceNumVerts(1));
         if obj.GaussPts.nNode == 1
            obj.N1 = obj.N1';
         end
      end

      function setQuad(obj,msh,GPoints)
         obj.mesh = msh;
         obj.GaussPts = GPoints;
         findLocDerBasisF(obj);
         findLocBasisF(obj);
         obj.detJ = zeros(1,obj.GaussPts.nNode);
      end
   end

end
