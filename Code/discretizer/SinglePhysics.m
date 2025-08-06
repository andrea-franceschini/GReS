classdef SinglePhysics < handle
   properties (Access = public)
      J
      rhs
   end

   properties
      model
      simParams
      dofm
      mesh
      elements
      faces
      material
      state
      fldId
      bcs
   end
   
   methods
      function obj = SinglePhysics(symmod,params,dofManager,grid,mat,bc,state)
         obj.model = symmod;
         obj.simParams = params;
         obj.dofm = dofManager;
         obj.mesh = grid.topology;
         obj.elements = grid.cells;
         obj.faces = grid.faces;
         obj.material = mat;
         obj.bcs = bc;
         obj.state = state;
         fld = feval([class(obj) '.getField']);
         obj.fldId = obj.dofm.getFieldId(fld);
      end

      function applyNeuBC(obj,dofs,vals)
         % Base Neumann BCs application method
         % Vals are subtracted since we solve du = J\(-rhs)!
         obj.rhs(dofs) = obj.rhs(dofs) - vals; 
      end

      function applyDirBC(obj,~,dofs,varargin)
         % Base application of boundary condition to jacobian block. This
         % method only implements standard Dirichlet BCs on diagonal
         % jacobian blocks, but can be overridden by specific
         % implementations.
         % This version works with incremental linear system (vals = 0).
         % BC values for Dirichlet BC are zero (since the
         % system is solved in incremental form)
         % set Dir rows to zero
         obj.J = obj.J';
         obj.J(:,dofs) = 0; % setting columns is much faster
         obj.J = obj.J';
         % Update rhs with columns to be removed
         %obj.rhs = obj.rhs - obj.J(:,dofs)*vals;
         % set obj.rhs vals to dir vals
         obj.rhs(dofs) = 0;
         obj.J(:,dofs) = 0;
         % modify diagonal entries of K (avoiding for loop and accessing)
         Jdiag = diag(obj.J);
         Jdiag(dofs) = 1;
         obj.J = obj.J - diag(diag(obj.J)) + diag(Jdiag);
      end

      function J = getJacobian(obj,varargin)
        J = obj.J;
      end

      function rhs = getRhs(obj,varargin)
        rhs = obj.rhs;
      end
   end

   methods (Static)
     function map = registerSolver(map,solverName)
       % update map of available solvers
       fld = eval([solverName '.getField']);
       map(fld) = str2func(solverName);
     end

   end
end