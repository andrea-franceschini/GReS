function [J,rhs] = applyBCFaults(NL,J,rhs,t)
% Apply Boundary condition to fault mechanics system
% Block dof indexing is used employing getContactDoF method of mortar
% faults class
% Impose BC to the linearized system (Jacobian matrix + RHS)
for domID = 1:NL.mortar.nDom
   model = NL.models(domID).ModelType;
   bound = NL.models(domID).BoundaryConditions;
   keys = bound.db.keys;
   for i = 1 : length(keys)
      dirVal = []; % if stays empty Penalty method is used
      rhsVal = []; % if stays empty Penalty method is used
      cond = bound.getCond(keys{i});
      type = bound.getType(keys{i});
      switch cond
         case 'NodeBC'
            bcDofs = bound.getDofs(keys{i});
            switch bound.getType(keys{i})
               case 'Neu'
                  rhsVal = - bound.getVals(keys{i}, t);
            end
         case 'SurfBC'
            if isFEMBased(model,'Poro')
               bcDofs = bound.getDofs(keys{i});
               switch type
                  case 'Neu'
                     entitiesInfl = bound.getEntitiesInfluence(keys{i});
                     q = bound.getVals(keys{i}, t);
                     rhsVal = - entitiesInfl*q;
               end
            end
      end

      % Map single domain dofs to global linear system
      % depending on the domain type
      if domID == NL.mortar.tagMaster
         dof = bcDofs;
      elseif domID == NL.mortar.tagSlave
         bcDofs = bcDofs + 3*NL.mortar.totNodMaster;
         if strcmp(type,'Dir') && ~strcmp(NL.multType,'P0') % remove constraint on slave dofs belonging to the interface
            dof = bcDofs(~ismember(bcDofs,NL.dofMap.intSlave));
         else % for P0 no need to remove constraint on slave dofs belonging to the interface
            dof = bcDofs;
         end
      end
      % mapping is trivial with only two domains
      switch type
         case 'Dir' % Dirichlet BC
            nrows = size(J,1);
            if isempty(rhsVal) && isempty(dirVal) % FEM Dirichlet BCs
               vals = zeros(numel(dof),1);
               [J,rhs] = applyDir(dof,vals,J,rhs);
            else
               J(nrows*(dof-1) + dof) = J(nrows*(dof-1) + dof) + dirVal;
               rhs(dof) = rhs(dof) + rhsVal;
            end
         case 'Neu'
            rhs(dof) = rhs(dof) + rhsVal;
         otherwise
            rhs(dof) = rhs(dof) + rhsVal;
      end
   end
end
end