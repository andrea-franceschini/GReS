classdef Triangle < handle
  % TRIANGLE element class

  properties (Access = private)
    mesh
%     vol
%     volSign
%     volNod
  end
  
  properties (Access = public)
    cellCentroid;
  end

  methods (Access = public)
    % Class constructor method
    function obj = Triangle(mesh)
      % Calling the function to set element data
      obj.setTriangle(mesh);
%       computeCellCentroid(obj);
    end
    %
    function [mat] = getDerBasisF(obj,el)
      %       for el = 1:obj.mesh.nCells
%       top = obj.mesh.cells(el,:);
      %         i = obj.mesh.cells(el,1);
      %         j = obj.mesh.cells(el,2);
      %         k = obj.mesh.cells(el,3);
      %         m = obj.mesh.cells(el,4);
      %
      inv_A = inv([1 obj.mesh.coordinates(obj.mesh.surfaces(el,1),1:2);
                   1 obj.mesh.coordinates(obj.mesh.surfaces(el,2),1:2);
                   1 obj.mesh.coordinates(obj.mesh.surfaces(el,3),1:2)]);
      mat = inv_A(2:3,:);
      %
      %         vol = (obj.mesh.coordinates(top,1)')*mat(2,:)';
      %       end
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
        vol(i) = det([1 obj.mesh.coordinates(top(1),:);
                      1 obj.mesh.coordinates(top(2),:);
                      1 obj.mesh.coordinates(top(3),:);])/2;
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

%     function v = getVolume(obj,el)
%       v = obj.vol(el);
%     end
%     
%     function vS = getVolumeSign(obj,el)
%       vS = obj.volSign(el);
%     end
  end
  
  methods (Access = private)
      function setTriangle(obj,mesh)
      obj.mesh = mesh;
    end
    
    function computeCellCentroid(obj)
      obj.cellCentroid = sparse(repelem(1:obj.mesh.nSurfaces,obj.mesh.surfNumVerts), ...
          nonzeros((obj.mesh.surfaces)'),repelem((obj.mesh.surfNumVerts).^(-1), ...
          obj.mesh.surfNumVerts),obj.mesh.nSurfaces,obj.mesh.nNodes) ...
          * obj.mesh.coordinates;
    end
  end
end