function [mat, rhs] = applyDir(dofs, vals, mat, rhs)
% Apply Dirichlet BCs to a linear system (no penalty)
% set Dir rows to zero
mat(dofs,:) = 0;
% Update rhs with columns to be removed
rhs = rhs - mat(:,dofs)*vals;
% set rhs vals to dir vals
rhs(dofs) = vals;
mat(:,dofs) = 0;
% modify diagonal entries of K (avoiding for loop and accessing)
Kdiag = diag(mat);
Kdiag(dofs) = 1;
mat = mat - diag(diag(mat)) + diag(Kdiag);
end