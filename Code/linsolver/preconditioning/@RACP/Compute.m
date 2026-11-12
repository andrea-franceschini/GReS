% Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
function A = Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);
   
   % Make the matrix 2x2
   A = obj.condenseDomains(A);

   % Compute the test space
   idxMain = obj.getKindBlocks("variable");
   TV0 = obj.buildTestSpace(obj.problemsolver, idxMain);

   % Compute and apply Ruiz diagonal scaling if requested
   A = obj.Ruiz.Compute(A);
   TV0 = obj.Ruiz.applyDinv(TV0,1);

   % Compute the augmented matrix
   [A11_aug,inv_D22] = obj.cpt_localAug(A{1,1},A{1,2},A{2,1},A{2,2},obj.PrecSym);
   
   % Compute the amg for block 11
   obj.AMG.Compute(A11_aug,obj.PrecSym,TV0,true);
  
   obj.Apply_L = @(x) obj.ApplyLeft(x,A11_aug,A{1,2},A{2,1},inv_D22);
   obj.Apply_R = @(x) obj.ApplyRight(x);
   
end

