 classdef Hexahedron < handle
  % HEXAHEDRON element class
 
  properties (Access = private)   % PRIVATE
%
% NODE ORDERING ASSUMPTION (same as Gmsh output):
% Hexahedron (8 nodes):             
% 
%        v
% 4----------3            
% |\     ^   |\           
% | \    |   | \          
% |  \   |   |  \         
% |   8------+---7        
% |   |  +-- |-- | -> u   
% 1---+---\--2   |        
%  \  |    \  \  |        
%   \ |     \  \ |        
%    \|      w  \|        
%     5----------6        

    coordLoc = [-1 -1 -1;
                 1 -1 -1;
                 1  1 -1;
                -1  1 -1;
                -1 -1  1;
                 1 -1  1;
                 1  1  1;
                -1  1  1]
    GaussPts
    J1
    mesh
    J
    detJ
    N1
  end

  methods (Access = public)
    % Class constructor method
    function obj = Hexahedron(msh,GPoints)
       obj.setHexahedron(msh,GPoints);
    end
       
    function [outVar1,outVar2] = getDerBasisFAndDet(obj,el,flOut)   % mat,dJWeighed
%       findJacAndDet(obj,el);  % OUTPUT: J and detJ
      % Find the Jacobian matrix of the isoparametric map and its determinant
      %
      % Possible ways of calling this function are:
      %    1) [mat,dJWeighed] = getDerBasisFAndDet(obj,el,1)
      %    2) mat = getDerBasisFAndDet(obj,el,2)
      %    3) dJWeighed = getDerBasisFAndDet(obj,el,3)
      obj.J = pagemtimes(obj.J1,obj.mesh.coordinates(obj.mesh.cells(el,:),:));
      if flOut == 3 || flOut == 1
%         obj.detJ = arrayfun(@(x) det(obj.J(:,:,x)),1:obj.GaussPts.nNode);
          for i=1:obj.GaussPts.nNode
            obj.detJ(i) = det(obj.J(:,:,i));
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
   
    
    function N1Mat = getBasisFinGPoints(obj)
      N1Mat = obj.N1;
    end
    
    function [vol,cellCentroid] = findVolumeAndCentroid(obj,idHexa)
      % Find the volume of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      vol = zeros(length(idHexa),1);
%       obj.volNod = zeros(obj.mesh.nNodes,1);
      cellCentroid = zeros(length(idHexa),3);
      i = 0;
      for el = idHexa
        i = i + 1;
        dJWeighed = getDerBasisFAndDet(obj,el,3);
        vol(i) = sum(dJWeighed);
        assert(vol(i)>0,'Volume less than 0');
        gPCoordinates = getGPointsLocation(obj,el);
        cellCentroid(i,:) = obj.detJ * gPCoordinates/vol(i);
      end
    end
    
    function nodeVol = findNodeVolume(obj,idHexa)
      nodeVol = zeros(8*length(idHexa),1);
      ptr = 0;
      for el = idHexa
        dJWeighed = obj.getDerBasisFAndDet(el,3);
        nodeVol(ptr+1:ptr+8) = obj.N1'*dJWeighed';
        ptr = ptr + 8;
      end
    end
    
    function gPCoordinates = getGPointsLocation(obj,el)
      % Get the location of the Gauss points in the element in the physical
      % space
      gPCoordinates = obj.N1*obj.mesh.coordinates(obj.mesh.cells(el,:),:);
    end
  end

  methods (Access = private)
    function findLocDerBasisF(obj)
      % Compute derivatives in the reference space for all Gauss points
      obj.J1 = zeros(obj.mesh.nDim,obj.mesh.cellNumVerts(1),obj.GaussPts.nNode);
      %
      % d(N)/d\csi
      d1 = bsxfun(@(i,j) 1/8*obj.coordLoc(j,1).* ...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)).* ...
          (1+obj.coordLoc(j,3).*obj.GaussPts.coord(i,3)), ...
          (1:obj.GaussPts.nNode)',1:obj.mesh.cellNumVerts(1));
      %
      % d(N)/d\eta
      d2 = bsxfun(@(i,j) 1/8*obj.coordLoc(j,2).* ...
          (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).* ...
          (1+obj.coordLoc(j,3).*obj.GaussPts.coord(i,3)), ...
          (1:obj.GaussPts.nNode)',1:obj.mesh.cellNumVerts(1));
      % d2 = 1/8.*coord_loc(:,2).*(1+coord_loc(:,1).*pti_G(1)).*(1+coord_loc(:,3).*pti_G(3));
      %
      % d(N)/d\zeta
      d3 = bsxfun(@(i,j) 1/8*obj.coordLoc(j,3).* ...
          (1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).* ...
          (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)), ...
          (1:obj.GaussPts.nNode)',1:obj.mesh.cellNumVerts(1));
      % d3 = 1/8.*coord_loc(:,3).*(1+coord_loc(:,1).*pti_G(1)).*(1+coord_loc(:,2).*pti_G(2));
      obj.J1(1,1:obj.mesh.cellNumVerts(1),1:obj.GaussPts.nNode) = d1';
      obj.J1(2,1:obj.mesh.cellNumVerts(1),1:obj.GaussPts.nNode) = d2';
      obj.J1(3,1:obj.mesh.cellNumVerts(1),1:obj.GaussPts.nNode) = d3';
    end
    
    function findLocBasisF(obj)
      % Find the value the basis functions take at the Gauss points
      obj.N1 = bsxfun(@(i,j) 1/8*(1+obj.coordLoc(j,1).*obj.GaussPts.coord(i,1)).* ...
                     (1+obj.coordLoc(j,2).*obj.GaussPts.coord(i,2)).* ...
                     (1+obj.coordLoc(j,3).*obj.GaussPts.coord(i,3)), ...
                     (1:obj.GaussPts.nNode)',1:obj.mesh.cellNumVerts(1));
    end
    
    function setHexahedron(obj,msh,GPoints)
      obj.mesh = msh;
      obj.GaussPts = GPoints;
      findLocDerBasisF(obj);
      findLocBasisF(obj);
      obj.detJ = zeros(1,obj.GaussPts.nNode);
    end
  end

  end
