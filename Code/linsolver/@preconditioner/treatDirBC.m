% Function for treating the dirichlet boundary conditions
function A = treatDirBC(obj,A,symMat)

   n11 = size(A{1,1},1);

   % Treat Dirichlet boundary conditions
   % Identify target indices
   D = sum(spones(A{1,1}));
   ind_dir_dof = find(D==1);
   ind_col_rem = sum(spones(A{1,2}))==1;
   [ind_dir_lag,~,~] = find(A{1,2}(:,ind_col_rem));
   ind_dir = union(ind_dir_dof,ind_dir_lag);

   % Native Column Zeroing 
   A{1,1}(:,ind_dir) = 0;
   A{2,1}(:,ind_dir) = 0;

   % A11 Row Zeroing
   A{1,1} = A{1,1}';
   A{1,1}(:,ind_dir) = 0;
   A{1,1} = A{1,1}';
   
   if symMat(1,2) == 1
      % If the matrix is symmetric then I can fix the dir copying the transposed block
      A{1,2} = A{2,1}';
   else
      % If the matrix is not symmetric I need to fix the 1,2 matrix itself
      A12t = A{1,2}';
      A12t(:,ind_dir) = 0;
      A{1,2} = A12t';
   end 

   % Diagonal Restoration
   fac = max(D);
   D_diag = zeros(n11,1);
   D_diag(ind_dir,1) = fac;
   A{1,1} = A{1,1} + spdiags(D_diag, 0, n11, n11);

end
