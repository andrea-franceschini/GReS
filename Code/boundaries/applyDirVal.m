function applyDirVal(bound, t, resState)
  % Apply Dirichlet conditions to the solution vector
  keys = bound.db.keys;
  for i = 1 : length(keys)
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
  end
end