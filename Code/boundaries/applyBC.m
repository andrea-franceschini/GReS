function applyBC(bound, t, syst)
  % Impose BC to the linearized system (Jacobian matrix + RHS)
  % The Penalty method is used for the Dirichlet conditions
  %
  maxVal = max(abs(syst.K), [], 'all');
  %
  keys = bound.db.keys;
  for i = 1 : length(keys)
    if strcmp(bound.getType(keys{i}), 'Neu')  % Apply Neumann conditions,if any
      syst.rhs(bound.getDofs(keys{i})) = syst.rhs(bound.getDofs(keys{i})) - bound.getVals(keys{i}, t);
    elseif strcmp(bound.getType(keys{i}), 'Dir')  % Apply Dirichlet conditions
      nrows = size(syst.K,1);
      syst.K(nrows*(bound.getDofs(keys{i})-1)+bound.getDofs(keys{i})) = maxVal*1.e10;
      syst.rhs(bound.getDofs(keys{i})) = 0;
    end
  end
end
