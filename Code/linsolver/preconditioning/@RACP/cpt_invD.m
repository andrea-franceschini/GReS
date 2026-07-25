function [inv_D22] = cpt_invD(obj,A11,A12,A21,A22)
   
   % Get the dimensions of the block
   n22 = size(A12,2);

   % Compute local augmentation
   aug = zeros(size(A22,1),1);
   D_11 = full(diag(A11));
   mean_diag_A = mean(D_11);

   A21_T = A21';
   for icol = 1:n22
      v12 = A12(:,icol);
      v21 = A21_T(:,icol);
      [ii_12,~,bb_12] = find(v12);
      [ii_21,~,bb_21] = find(v21);
      if (numel(ii_12)+numel(ii_21) > 0)
         BB = bb_12*bb_21';
         AA = A11(ii_12,ii_21);
         m_a = max(eig(full(AA)));
         m_b = max(eig(full(BB)));
         
         if m_a == 0
            m_a = mean_diag_A;
         end
         alpha = m_a / m_b;
         aug(icol) = 1 / alpha;
      end
   end

   aug_mat = diag(sparse(aug));
   A22_aug = A22 - obj.gamma*aug_mat;

   % Compute the inverse of the 22 augmentation
   inv_D22 = -inv(diag(diag(A22_aug)));
end