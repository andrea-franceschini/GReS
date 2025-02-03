classdef CouplingPhysics < handle
   % General coupling solver   
   properties 
      fields
      J = cell(2,1)             % 2x1 cell with jacobian blocks of fields
      rhs = cell(2,1)           % 2x1 cell array with rhs blocks of fields
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
      function obj = CouplingPhysics(field1,field2,symmod,params,dofManager,grid,mat,data)
         obj.fields = {field1,field2};
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
      
      function applyDirBC(obj,field,dofs,varargin)
         % standard BC application for coupling blocks. 
         % This method implements only the Dirichlet BC application
         % It is ovverridden by subclasses: e.g Biot does not call this
         % method if flow is discretized with TPFA 
         Jid = strcmp(obj.fields,field);
         % set obj.rhs vals to dir vals
         obj.rhs{Jid}(dofs) = 0;
         % set row of J to 0
         obj.J{Jid} = (obj.J{Jid})';
         obj.J{Jid}(:,dofs) = 0;
         obj.J{Jid} = (obj.J{Jid})';
         % set column to 0
         obj.J{~Jid}(:,dofs) = 0;     
      end

      function applyNeuBC(varargin)
         % Neumann BCs are not applied to Coupling solvers
         return
      end

      function J = getJacobian(obj,fld)
         id = find(strcmp(obj.fields,fld));
         assert(any(id),['Invalid field identifier %s for Coupling Solver' ...
            'class'],fld);
         J = obj.J{id};
      end

      function rhs = getRhs(obj,fld)
         id = strcmp(obj.fields,fld);
         assert(any(id),['Invalid field identifier %s for Coupling Solver' ...
            'class'],fld);
         rhs = obj.rhs{id};
      end
   end
end

