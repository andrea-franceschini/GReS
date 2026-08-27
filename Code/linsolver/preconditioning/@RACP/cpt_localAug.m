function [A11_aug,inv_D22] = cpt_localAug(obj,A11,A12,A21,A22,symm)   
   simple_flag = false;

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
         if simple_flag
            m_a= max(D_11(ii_12)); 
            m_b = max(diag(BB));
         else
            AA = A11(ii_12,ii_21);
            m_a = max(eig(full(AA)));
            m_b = max(eig(full(BB)));
         end
         if m_a == 0
            m_a = mean_diag_A;
         end
         alpha = m_a / m_b;
         aug(icol) = 1 / alpha;
      end
   end

   aug_mat = diag(sparse(aug));
   A22_aug = A22 - obj.gamma*aug_mat;

   % Compute augmented 11 block
   inv_D22 = -inv(diag(diag(A22_aug)));
   ADD = A12*inv_D22*A21; 

   % Strong Symmetrization if the matrix was seen as symmetric
   if symm == 1
      ADD = 0.5*(ADD+ADD');
   end
   A11_aug = A11+ADD;

end