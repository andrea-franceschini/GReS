classdef DofManagerSedMatrix < handle
  
  properties
    dof (:,:,:) uint64
    ndofs uint64
  end

  properties(Access = private)
    dim uint16
  end

  methods
    function obj = DofManagerSedMatrix(ncells)
      obj.dim = ncells;
      obj.ndofs = prod(ncells);
      obj.dof = reshape(1:obj.ndofs,ncells);
    end

    function dof = getFromIJK(obj,i,j,k)
      dof = obj.dof(i,j,k);
    end

    function neigh = findNeighbors(obj)
      neigh = zeros(obj.ndofs,6,'uint64');
      for k = 1:obj.dim(3)
        for i = 1:obj.dim(1)
          for j = 1:obj.dim(2)
            cellId = obj.dof(i, j, k);
            % Neighbor Mapping:
            if i > 1
              neigh(cellId, 3) = obj.dof(i-1,j,k);
            end
            if i < obj.dim(1)
              neigh(cellId, 4) = obj.dof(i+1,j,k);
            end
            if j > 1
              neigh(cellId, 1) = obj.dof(i,j-1,k);
            end
            if j < obj.dim(2)
              neigh(cellId, 2) = obj.dof(i,j+1,k);
            end
            if k > 1
              neigh(cellId, 5) = obj.dof(i,j,k-1);
            end
            if k < obj.dim(3)
              neigh(cellId, 6) = obj.dof(i,j,k+1);
            end
          end
        end
      end
    end

    function dof = getDof(obj)
      dof=obj.dof(:);
    end



  end
end