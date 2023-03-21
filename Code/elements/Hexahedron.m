classdef Hexahedron < handle
  % HEXAHEDRON element class

  properties (Access = public)
    vol
    volNod
  end
  
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
%     vol
    J
    detJ
%     invJ;
    N1
  end

  methods (Access = public)
    % Class constructor method
    function obj = Hexahedron(msh,GPoints)
       % Calling the function to set element data
%        obj.setElementData(data);
       % ALLOCATE j
       obj.setHexahedron(msh,GPoints);
       findVolume(obj);
%        obj.J = zeros(3,3,obj.GaussPts.nNode);
%        obj.detJ = zeros(obj.GaussPts.nNode,1);
    end
       
    function [mat,dJWeighed] = getDerBasisFAndDet(obj,el)
      findJacAndDet(obj,el);  % OUTPUT: J and detJ
      invJTmp = arrayfun(@(x) inv(obj.J(:,:,x)),1:obj.GaussPts.nNode,'UniformOutput',false);
      obj.J = reshape(cell2mat(invJTmp),obj.mesh.nDim,obj.mesh.nDim,obj.GaussPts.nNode); %possibly we can overwrite J
      clear invJTmp
      mat = pagemtimes(obj.J,obj.J1);
      dJWeighed = obj.detJ.*(obj.GaussPts.weight)';
    end
    
    function findJacAndDet(obj,el)
      % Find the Jacobian matrix of the isoparametric map and its determinant
      obj.J = pagemtimes(obj.J1,obj.mesh.coordinates(obj.mesh.cells(el,:),:));
      obj.detJ = arrayfun(@(x) det(obj.J(:,:,x)),1:obj.GaussPts.nNode);
    end
    
    function coord = getGPointsLocation(obj,el)
      % Get the location of the Gauss points in the element in the physical
      % space
      coord = obj.N1*obj.mesh.coordinates(obj.mesh.cells(el,:),:);
    end
    
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
    
    function v = getVolume(obj,el)
      v = obj.vol(el);
    end
  end

  methods (Access = private)
    function findVolume(obj)
      % Find the volume of the cells using the determinant of the Jacobian
      % of the isoparameric transformation
      obj.vol = zeros(obj.mesh.nCells,1);
      obj.volNod = zeros(obj.mesh.nNodes,1);
      for el=1:obj.mesh.nCells
        findJacAndDet(obj,el)
        obj.vol(el) = (obj.detJ) * (obj.GaussPts.weight);
        assert(obj.vol(el)>0,'Volume less than 0');
        top = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
        for i=1:obj.mesh.cellNumVerts(el)
          obj.volNod(top(i)) = obj.volNod(top(i)) + obj.vol(el)/obj.mesh.cellNumVerts(el);
        end
      end
    end
    
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
    end
  end

end