function applyDirVal(mod, bound, t, resState)
  % Apply Dirichlet conditions to the solution vector
  keys = bound.db.keys;
  for i = 1 : length(keys)
    %Nodal dirichlet conditions
    if (strcmp(bound.getCond(keys{i}),"NodeBC") && ...
       strcmp(bound.getType(keys{i}), 'Dir')) || ...
       strcmp(bound.getCond(keys{i}), 'ElementBC')
%       if strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
      if strcmp(bound.getPhysics(keys{i}), 'Flow')
        resState.pressure(bound.getEntities(keys{i})) = bound.getVals(keys{i}, t);
      elseif strcmp(bound.getPhysics(keys{i}), 'Poro')
        resState.dispConv(bound.getEntities(keys{i})) = bound.getVals(keys{i}, t);
        resState.dispCurr(bound.getEntities(keys{i})) = bound.getVals(keys{i}, t);
      end
%       end
    end

    %surface Dirichlet conditions
    if (strcmp(bound.getCond(keys{i}),"SurfBC") && ...
       strcmp(bound.getType(keys{i}), 'Dir')) 
       if strcmp(bound.getPhysics(keys{i}), 'Flow') && isFEMBased(mod,'Flow')
          surfVal = bound.getVals(keys{i}, t);
          nodMap = bound.getEntitiesInfluence(keys{i});
          resState.pressure(bound.getCompLoadedEntities(keys{i})) = nodMap*surfVal;
       elseif strcmp(bound.getPhysics(keys{i}), 'Poro') && isFEMBased(mod,'Poro')
          surfVal = bound.getVals(keys{i}, t);
          nodMap = bound.getEntitiesInfluence(keys{i});
          resState.dispConv(bound.getCompLoadedEntities(keys{i})) = nodMap*surfVal;
          resState.dispCurr(bound.getCompLoadedEntities(keys{i})) = nodMap*surfVal;
       end
    end
  end