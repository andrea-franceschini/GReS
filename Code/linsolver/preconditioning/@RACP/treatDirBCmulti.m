% Function for treating Dirichlet boundary conditions in an NxN block system
function A = treatDirBCmulti(obj, A, sym)
   N = size(A, 1); % Dimension of the block matrix
   
   % Loop through primal blocks (1 to N-1), ignoring the N,N block
   for i = 1:N-1
       n_ii = size(A{i,i}, 1);
       
       % Identify target indices for block i
       D = sum(spones(A{i,i}), 1);
       ind_dir_dof = find(D == 1);
       
       % Lagrange multiplier coupling with constraint block N
       if ~isempty(A{i,N})
           ind_col_rem = sum(spones(A{i,N}), 1) == 1;
           [ind_dir_lag, ~, ~] = find(A{i,N}(:, ind_col_rem));
           ind_dir = union(ind_dir_dof, ind_dir_lag);
       else
           ind_dir = ind_dir_dof;
       end
       
       % Skip if no Dirichlet DOFs are present for this block
       if isempty(ind_dir)
           continue;
       end
       
       % 1. Native Column Zeroing across all row blocks j for column block i
       for j = 1:N
           if ~isempty(A{j,i})
               A{j,i}(:, ind_dir) = 0;
           end
       end
       
       % 2. Diagonal Block A_{i,i} Row Zeroing
       A{i,i} = A{i,i}';
       A{i,i}(:, ind_dir) = 0;
       A{i,i} = A{i,i}';
       
       % 3. Diagonal Restoration for A_{i,i}
       fac = max(D);
       if isempty(fac) || fac == 0
           fac = 1; % Fallback scale factor
       end
       D_diag = zeros(n_ii, 1);
       D_diag(ind_dir, 1) = fac;
       A{i,i} = A{i,i} + spdiags(D_diag, 0, n_ii, n_ii);
       
       % 4. Off-Diagonal Row Zeroing / Symmetry Enforcement for A_{i,j} (j ~= i)
       for j = 1:N
           if j == i
               continue;
           end
           
           % Check symmetry condition: A_{i,i} is symmetric, A_{i,j} = A_{j,i}', and A_{j,j} is symmetric
           % (When j = N, this checks A_{i,i}, A_{i,N}/A_{N,i}, and the A_{N,N} block)
           is_symmetric = sym(i,i) && sym(i,j) && sym(j,j);
           
           if is_symmetric
               % Copy the transposed block to enforce symmetry directly
               A{i,j} = A{j,i}';
           else
               % Zero out rows of A_{i,j} corresponding to ind_dir
               if ~isempty(A{i,j})
                   A_ij_t = A{i,j}';
                   A_ij_t(:, ind_dir) = 0;
                   A{i,j} = A_ij_t';
               end
           end
       end
   end
end