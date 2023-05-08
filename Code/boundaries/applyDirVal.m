function applyDirVal(bound, t, resState)
  % Apply Dirichlet conditions to the solution vector
  keys = bound.db.keys;
  for i = 1 : length(keys)
    if (ismember(bound.getCond(keys{i}), ["NodeBC","SurfBC"]) && ...
       strcmp(bound.getType(keys{i}), 'Dir')) || ...
       strcmp(bound.getCond(keys{i}), 'ElementBC')
%       if strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
      if strcmp(bound.getPhysics(keys{i}), 'Flow')
        resState.pressure(bound.getDofs(keys{i})) = bound.getVals(keys{i}, t);
      elseif strcmp(bound.getPhysics(keys{i}), 'Poro')
        resState.displ(bound.getDofs(keys{i})) = bound.getVals(keys{i}, t);
      end
%       end
    end
  end
end