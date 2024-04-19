function [K, f] = applyDir(dofs, vals, K, f)
% Apply Dirichlet BCs with an improved algorithm compared to Penalty
% approach
% set Dir rows to zero
K(dofs,:) = 0;
% Update rhs with columns to be removed
f = f - K(:,dofs)*vals;
% set rhs vals to dir vals
f(dofs) = vals;
K(:,dofs) = 0;
for i = 1:length(dofs)
    K(dofs(i),dofs(i)) = 1;
end

