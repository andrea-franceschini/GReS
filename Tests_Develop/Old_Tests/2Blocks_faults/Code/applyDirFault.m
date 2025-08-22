function applyDirFault(NL)
slaveNodesID = NL.mortar.idSlave;
% Apply Dirichlet conditions to the solution vector in model with faults
for domID = 1:NL.mortar.nDom
   bound = NL.models(domID).BoundaryConditions;
   keys = bound.db.keys;
   for i = 1:length(keys)
      if ~strcmp(bound.getType(keys{i}), 'Dir') || ~(strcmp(bound.getPhysics(keys{i}), 'Poro'))
         continue % discard non dirichlet BCs or non poro conditions
      end
      % get entities and BC dofs
      if strcmp(bound.getCond(keys{i}),"NodeBC")
         dofBC = bound.getEntities(keys{i});
         % component multiplication
         dofBC = bound.getCompEntities(keys{i},dofBC);
         valBC = bound.getVals(keys{i}, NL.t);
      elseif strcmp(bound.getCond(keys{i}),"SurfBC")
         dofBC = bound.getCompLoadedEntities(keys{i});
         surfVal = bound.getVals(keys{i}, NL.t);
         nodMap = bound.getEntitiesInfluence(keys{i});
         valBC = nodMap*surfVal;
      end
      % remove interface slave dofs from BC dofs
      if domID == NL.mortar.tagSlave
         id = ~ismember(dofBC,get_dof(slaveNodesID));
         dofBC = dofBC(~ismember(dofBC,get_dof(slaveNodesID)));
         valBC = valBC(id);
      end
      NL.state(domID).dispConv(dofBC) = valBC;
      NL.state(domID).dispCurr(dofBC) = valBC;
   end
end
end