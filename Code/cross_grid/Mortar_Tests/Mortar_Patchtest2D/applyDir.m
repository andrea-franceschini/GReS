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
% modify diagonal entries of K (avoiding for loop and accessing)
Kdiag = diag(K);
Kdiag(dofs) = 1;
K = K - diag(diag(K)) + diag(Kdiag);
end

