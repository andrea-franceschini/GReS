classdef SinglePhysicSolver < handle
   properties
      J
      rhs
      field
      model
      simParams
      dofm
      mesh
      elements
      faces
      material
      GaussPts
   end
   
   methods
      function obj = SinglePhysicSolver(field,symmod,params,dofManager,grid,mat,data)
         obj.field = field;
         obj.model = symmod;
         obj.simParams = params;
         obj.dofm = dofManager;
         obj.mesh = grid.topology;
         obj.elements = grid.cells;
         obj.faces = grid.faces;
         obj.material = mat;
         if ~isempty(data)
            obj.GaussPts = data{1};
         end
      end
      
      function applyBC(obj)
         % Base application of boundary condition to jacobian block This
         % method only implements standard Dirichlet BCs on diagonal
         % jacobian block but can be overridden by specific implementations
      end
   end
end

