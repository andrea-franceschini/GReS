function [J, rhs] = applyDirBlock(locDofs,vals,J,rhs,id)
% Apply Dirichlet BCs to a block linear system
n = size(J,1);
% set dir rows to 0
for i = 1:n
    J{id,i}(locDofs,:)=0;
end
% Update rhs with columns to be removed
for i = 1:n
    rhs{i} = rhs{i} - J{i,id}(:,locDofs)*vals;
end
% set rhs vals to dir vals
rhs{id}(locDofs) = vals;
% set Dir columns to 0
for i = 1:n
    J{i,id}(:,locDofs)=0;
end
% modify diagonal entries of J
Kdiag = diag(J{id,id});
Kdiag(locDofs) = 1;
J{id,id} = J{id,id} - diag(diag(J{id,id})) + diag(Kdiag);
end