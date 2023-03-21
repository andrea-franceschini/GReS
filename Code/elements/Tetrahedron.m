classdef Tetrahedron < handle
  % TETRAHEDRON element class

  properties (Access = public)
    mesh
    vol
    volSign
    volNod
  end
  
  properties (Access = public)
    cellCentroid;
  end

  methods (Access = public)
    % Class constructor method
    function obj = Tetrahedron(mesh)
      % Calling the function to set element data
      obj.setTetrahedron(mesh);
      computeCellCentroid(obj);
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
      inv_A = inv([1 obj.mesh.coordinates(obj.mesh.cells(el,1),:);
                   1 obj.mesh.coordinates(obj.mesh.cells(el,2),:);
                   1 obj.mesh.coordinates(obj.mesh.cells(el,3),:);
                   1 obj.mesh.coordinates(obj.mesh.cells(el,4),:)]);
      mat = inv_A(2:4,:);
      %
      %         vol = (obj.mesh.coordinates(top,1)')*mat(2,:)';
      %       end
    end

    function v = getVolume(obj,el)
      v = obj.vol(el);
    end
    
    function vS = getVolumeSign(obj,el)
      vS = obj.volSign(el);
    end
  end
   
  methods (Access = private)
%   Elements volume calculation
    function findVolume(obj)
      obj.vol = zeros(obj.mesh.nCells,1);
      obj.volSign = ones(obj.mesh.nCells,1);
      obj.volNod = zeros(obj.mesh.nNodes,1);
      for el = 1:obj.mesh.nCells
        top = obj.mesh.cells(el,1:obj.mesh.cellNumVerts(el));
        obj.vol(el) = det([1 obj.mesh.coordinates(top(1),:);
                           1 obj.mesh.coordinates(top(2),:);
                           1 obj.mesh.coordinates(top(3),:);
                           1 obj.mesh.coordinates(top(4),:)])/6;
        if obj.vol(el) < 0
          obj.volSign(el) = -1;
          obj.vol(el) = -obj.vol(el);
        end
        for i=1:obj.mesh.cellNumVerts(el)
          obj.volNod(top(i)) = obj.volNod(top(i)) + obj.vol(el)/obj.mesh.cellNumVerts(el);
        end
      end
    end
  end
  
  methods (Access = private)
    function setTetrahedron(obj,mesh)
      obj.mesh = mesh;
      if any(obj.mesh.cellNumVerts == 10)
        error('Quadratic tetrahedron not available');
      end
      findVolume(obj);
    end
    
    function computeCellCentroid(obj)
      obj.cellCentroid = sparse(repelem(1:obj.mesh.nCells,obj.mesh.cellNumVerts), ...
          nonzeros((obj.mesh.cells)'),repelem((obj.mesh.cellNumVerts).^(-1), ...
          obj.mesh.cellNumVerts),obj.mesh.nCells,obj.mesh.nNodes) ...
          * obj.mesh.coordinates;
    end
  end
end