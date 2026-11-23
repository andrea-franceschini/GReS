% Function to compute the RACP preconditioner for the lagrange multiplier case (single physics multi domain)
function A = Compute(obj,A,symMat,varargin)

   if ~iscell(A)
      error("Parent preconditioner not configured to run with multiphys multidomain");
   end

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);
   
   % Compute the symmetry of the augmented matrix
   symMatAug = computeSymMatAug(obj,A,symMat);

   % Retrieve or compute cached unscaled test space
   if nargin >= 4 && ~isempty(varargin{1})
      TV0 = varargin{1};
   else
      TV0 = obj.getBlockTestSpace(obj.problemsolver);
   end

   % Compute and apply Ruiz diagonal scaling on full block matrix A if requested
   A = obj.Ruiz.Compute(A);
   TV0 = obj.Ruiz.applyDinv(TV0);

   % Condense only saddle-point interfaces into block K+1, keeping variable blocks separate
   A = obj.condenseDomains(A);

   % Compute augmented variable blocks directly at cell block level
   [A11_aug,inv_D22,A12,A21] = obj.cpt_localAug(A,obj.PrecSym);
   
   % Slice test space for variable blocks
   idxMain = obj.getKindBlocks("variable");
   TV0_main = TV0(idxMain);

   % Compute the inner preconditioner for augmented variable blocks
   obj.innerPrec.Compute(A11_aug,symMatAug(idxMain,idxMain),TV0_main,true);
  
   obj.Apply_L = @(x) obj.ApplyLeft(x,A12,A21,inv_D22);
   obj.Apply_R = @(x) obj.ApplyRight(x);
   
end




% Compute the symmetry of the augmented matrix correctly
function symMatAug = computeSymMatAug(obj,A,symMat)
   
   % Get the blocks to loop on
   idxMain = obj.getKindBlocks("variable");
   idxSupp = obj.getKindBlocks("interface");

   % Initialize the augmented matrix symmetry
   symMatAug = symMat;

   % Loop over the main blocks
   for i = idxMain
      % Loop over the associated interface blocks
      for j = idxSupp
         % If the interface matrix is nonzero and the main block is still
         % symmetric check if the augmentation will turn it nonsymmetric
         if nnz(A{i,j}) > 0 && symMatAug(i,i) ~= 0
            symMatAug(i,i) = floor(sum(symMat([i j],[i j]),'all')/4);
         end
      end
   end
end