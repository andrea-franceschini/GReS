% Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
function Compute(obj,A,sym,varargin)

   obj.PrecSym = sym;
   simple_flag = false;
   
   % Treat the multiple domains as if they were one and then use RACP
   if numel(A) ~= 4
      % checks if any of the matrix entries are 0x0 blocks
      [ZeroSpRow,ZeroSpCol] = find(cellfun(@(x) isempty(x), A));
      if ~isempty(ZeroSpRow)
         % get the correct number of rows
         rows = zeros(size(A,1),1);
         for i = 1:size(A,1)
            rows(i) = size(A{i,1},1);
         end 
         % assign the correct dimension to the matrices
         for i = 1:length(ZeroSpRow)
            A{ZeroSpRow(i),ZeroSpCol(i)} = sparse(rows(ZeroSpRow(i)),rows(ZeroSpCol(i)));
         end
      end
      nn = size(A,1);
      idxMain = 1:nn-obj.nInt;
      idxSupp = nn+1-obj.nInt:nn;
      A11 = cell2matrix(A(idxMain,idxMain));
      A12 = cell2matrix(A(idxMain,idxSupp));
      A21 = cell2matrix(A(idxSupp,idxMain));
      A22 = cell2matrix(A(idxSupp,idxSupp)); 
      
      clear A;
  
      A = {A11, A12; A21 A22};
   end

   % Get the dimensions of the block
   n22 = size(A{1,2},2);

   % If block 22 has dim 0 then resize it
   if size(A{2,2},1) ~= n22
      A{2,2} = sparse(n22,n22);
   end

   A = obj.treatDirBC(A,sym);

   % Compute local augmentation
   aug = zeros(size(A{2,2},1),1);
   D_11 = full(diag(A{1,1}));
   mean_diag_A = mean(D_11);

   A21_T = A{2,1}';
   for icol = 1:n22
      v12 = A{1,2}(:,icol);
      v21 = A21_T(:,icol);
      [ii_12,~,bb_12] = find(v12);
      [ii_21,~,bb_21] = find(v21);
      if (numel(ii_12)+numel(ii_21) > 0)
         BB = bb_12*bb_21';
         if simple_flag
            m_a= max(D_11(ii_12)); 
            m_b = max(diag(BB));
         else
            AA = A{1,1}(ii_12,ii_21);
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
   A22_aug = A{2,2} - obj.gamma*aug_mat;

   % Compute augmented 11 block
   inv_D22 = -inv(diag(diag(A22_aug)));
   ADD = A{1,2}*inv_D22*A{2,1}; 

   % Strong Symmetrization if the matrix was seen as symmetric
   if sym == 1
      ADD = 0.5*(ADD+ADD');
   end
   A11_aug = A{1,1}+ADD;
   
   % Compute the test space
   if(obj.phys == 0) % fluids
      TV0 = [];
      for i = 1:obj.nDom - obj.nInt
         TV = ones(size(A{i,1},1),1);
         TV0 = [TV0;TV];
      end
   elseif(obj.phys == 1 || obj.phys == 1.1) % true contact mechanichs physics is 1.1, general poromechanics is 1
      TV0 = [];
      for i = 1:obj.nDom
         TV = mk_rbm_3d(obj.domain(i).grid.coordinates);
         TV0 = [TV0;TV];
      end
      % if obj.DEBUGflag
      %    obj.TV0 = TV0;
      % end
   end

   % Compute the amg for block 11
   obj.AMG.Compute(A11_aug,sym,TV0,true);

   obj.Apply_L = @(x) obj.ApplyLeft(x,A11_aug,A{1,2},A{2,1},inv_D22);
   obj.Apply_R = @(x) obj.ApplyRight(x);
end

