function [globalsymm,maxval,symMat] = checkSymmetry(A,eps1)
% CHECKSYMMETRY Checks symmetry of numeric or block-cell matrix A against tolerance eps1

   % Base Case: Numeric Matrix
   if ~iscell(A)

      diffnorm = norm(A-A','f');
      Anorm = norm(A,'f');
      
      if diffnorm == 0 || Anorm == 0
         globalsymm = 1;
         maxval = 0;
      else
         relNorm = diffnorm/Anorm;

         if relNorm < eps1
            maxval = 0;
            globalsymm = 1;
         else
            maxval = relNorm;
            globalsymm = 0;
         end
      end
      symMat = globalsymm;
      return
   end
   
   % Allocate the stuff
   N = size(A,1);
   symm = ones(sum(1:N),1);
   val = zeros(sum(1:N),1);
   cont = 1;
   
   % Loop over the blocks
   for j = 1:N
      for i = 1:j
         if i == j
            % Diagonal Block
            [symm(cont), val(cont)] = checkSymmetry(A{i,i},eps1);
         elseif nnz(A{i,j}) > 0 || nnz(A{j,i}) > 0
            % Off-Diagonal Block
            diffnorm = norm(A{i,j}-A{j,i}','f');
            Anorm = 0.5 * (norm(A{i,j},'f') + norm(A{j,i},'f'));
            
            if diffnorm == 0 || Anorm == 0
               symm(cont) = 1;
               val(cont) = 0;
            else
               relNorm = diffnorm/Anorm;

               if relNorm < eps1
                  symm(cont) = 1;
                  val(cont) = 0;
               else
                  symm(cont) = 0;
                  val(cont) = relNorm;
               end
            end
         end
         cont = cont + 1;
      end
   end
   
   % Global mapping and output
   symMat = zeros(N, N);
   symMat(triu(true(N))) = symm;
   symMat = symMat + triu(symMat, 1).';
   
   globalsymm = min(symm);
   maxval = max(val);
end
