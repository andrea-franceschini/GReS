function state = applyDirVal(mod, bound, t, state)
  % Apply Dirichlet conditions to the solution vector
  keys = bound.db.keys;
  for i = 1 : length(keys)
    %Nodal dirichlet conditions
    if (strcmp(bound.getCond(keys{i}),"NodeBC") && ...
       strcmp(bound.getType(keys{i}), 'Dir')) || ...
       strcmp(bound.getCond(keys{i}), 'ElementBC')
%       if strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
      if strcmp(bound.getPhysics(keys{i}), 'Flow')
        state.pressure(bound.getEntities(keys{i})) = bound.getVals(keys{i}, t);
      elseif strcmp(bound.getPhysics(keys{i}), 'Poro')
        state.dispConv(bound.getEntities(keys{i})) = bound.getVals(keys{i}, t);
        state.dispCurr(bound.getEntities(keys{i})) = bound.getVals(keys{i}, t);
      end
%       end
    end

    %surface Dirichlet conditions
    if (strcmp(bound.getCond(keys{i}),"SurfBC") && ...
       strcmp(bound.getType(keys{i}), 'Dir')) 
       if strcmp(bound.getPhysics(keys{i}), 'Flow') && isFEMBased(mod,'Flow')
          surfVal = bound.getVals(keys{i}, t);
          nodMap = bound.getEntitiesInfluence(keys{i});
          state.pressure(bound.getCompLoadedEntities(keys{i})) = nodMap*surfVal;
       elseif strcmp(bound.getPhysics(keys{i}), 'Poro') && isFEMBased(mod,'Poro')
          surfVal = bound.getVals(keys{i}, t);
          nodMap = bound.getEntitiesInfluence(keys{i});
          state.dispConv(bound.getCompLoadedEntities(keys{i})) = nodMap*surfVal;
          state.dispCurr(bound.getCompLoadedEntities(keys{i})) = nodMap*surfVal;
       end
    end
  end