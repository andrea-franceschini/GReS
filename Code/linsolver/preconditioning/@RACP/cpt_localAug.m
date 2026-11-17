function [A11_aug,inv_D22,A12,A21] = cpt_localAug(obj,A,symm,varargin)
   simple_flag = false;

   if nargin >= 5
      % Legacy signature: cpt_localAug(obj, A11, A12, A21, A22, symm)
      A11 = A;
      A12 = symm;
      A21 = varargin{1};
      A22 = varargin{2};
      if numel(varargin) >= 3
         symm = varargin{3};
      else
         symm = 1;
      end
      isCellBlock = false;
   else
      % Block signature: cpt_localAug(obj, A, symm) where A is (K+1) x (K+1) cell array
      K = size(A, 1) - 1;
      A12 = cell2matrix(A(1:K, K+1));
      A21 = cell2matrix(A(K+1, 1:K));
      A22 = A{K+1, K+1};
      A11 = cell2matrix(A(1:K, 1:K));
      isCellBlock = true;
   end

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
   
   iscomplex = ~isreal(aug);

   if any(iscomplex)
     compAug = aug(iscomplex);
     r = real(compAug);
     c = imag(compAug);

     n = norm(c./r,'inf');

     if n > 100*eps
       error('Coumplex found in local augmentation')
     else
       aug = real(aug);
     end
   end

   aug_mat = diag(sparse(aug));
   A22_aug = A22 - obj.gamma*aug_mat;

   % Compute augmented 11 block
   inv_D22 = -inv(diag(diag(A22_aug)));

   if isCellBlock
      % Compute augmented blocks directly at the cell block level
      A11_aug = A(1:K, 1:K);
      for i = 1:K
         for j = 1:K
            if nnz(A{i, K+1}) > 0 && nnz(A{K+1, j}) > 0
               ADD_ij = A{i, K+1} * inv_D22 * A{K+1, j};
               if symm == 1 && i == j
                  ADD_ij = 0.5 * (ADD_ij + ADD_ij');
               end
               A11_aug{i, j} = A11_aug{i, j} + ADD_ij;
            end
         end
      end
      if symm == 1
         for i = 1:K
            for j = (i+1):K
               if nnz(A{i, K+1}) > 0 && nnz(A{K+1, j}) > 0
                  A11_aug{i, j} = 0.5 * (A11_aug{i, j} + A11_aug{j, i}');
                  A11_aug{j, i} = A11_aug{i, j}';
               end
            end
         end
      end
   else
      ADD = A12*inv_D22*A21; 
      if symm == 1
         ADD = 0.5*(ADD+ADD');
      end
      A11_aug = A11+ADD;
   end

end