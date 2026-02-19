% Function for treating the dirichlet boundary conditions
function treatDirBC(obj,A)
   % Get the dimensions of the blocks
   n11 = size(A{1,1},1);

   % Treat Dirichlet boundary conditions
   A{1,1} = A{1,1}';
   D = sum(spones(A{1,1}));
   ind_dir_dof = find(D==1);
   ind_col_rem = find(sum(spones(A{1,2}))==1);
   [ind_dir_lag,~,~] = find(A{1,2}(:,ind_col_rem));
   ind_dir = union(ind_dir_dof,ind_dir_lag);
   A{1,1}(:,ind_dir) = 0;
   A{1,1} = A{1,1}';
   A{2,1}(:,ind_dir) = 0;
   A{1,2} = A{2,1}';
   fac = max(D);
   D = zeros(n11,1);
   D(ind_dir,1) = fac;
   A{1,1} = A{1,1} + diag(sparse(D));
end
