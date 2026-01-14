% Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
function computeRACP(obj,A)

   simple_flag = false;

   % Get the dimensions of the block
   n22 = size(A{1,2},2);

   % If block 22 has dim 0 then resize it
   if size(A{2,2},1) ~= n22
      A{2,2} = sparse(n22,n22);
   end

   % Treat Dirichlet boundary conditions
   obj.treatDirBC(A);

   % Set RACP Gamma to 1
   gamma = 1.0;

   % Compute local augmentation
   aug = zeros(size(A{2,2},1),1);
   D_11 = full(diag(A{1,1}));
   mean_diag_A = mean(D_11);
   D_22 = full(diag(A{2,2}));
   A21_scaled_T = A{2,1}';
   for icol = 1:n22
      v12 = A{1,2}(:,icol);
      v21 = A21_scaled_T(:,icol);
      [ii_12,~,bb_12] = find(v12);
      [ii_21,~,bb_21] = find(v21);
      if (numel(ii_12)+numel(ii_21) > 0);
         BB = bb_12*bb_21';
         if simple_flag
            m_a= max(D_11(ii_12));
            m_b = max(diag(BB));
         else
            AA = A{1,1}(ii_12,ii_21);
            m_a = max(eig(full(AA)));
            m_b = max(eig(full(BB)));
         end
         m_a_sav = m_a;
         if m_a == 0
            m_a = mean_diag_A;
         end
         alpha = m_a / m_b;
         aug(icol) = 1 / alpha;
      end
   end

   aug_mat = diag(sparse(aug));
   A22_aug = A{2,2} - gamma*aug_mat;

   % Compute augmented 11 block
   inv_D22 = -inv(diag(diag(A22_aug)));
   ADD = A{1,2}*inv_D22*A{2,1}; ADD = 0.5*(ADD+ADD');
   A11_aug = A{1,1}+ADD;
   
   % For now impose the amg
   obj.PrecType = 'amg';

   % Compute the amg for block 11
   obj.computePrec(A11_aug);

   obj.MfunL = @(x) apply_RevAug(obj.Prec,A11_aug,A{1,2},inv_D22,x);
   obj.MfunR = @(x) x;
end
