% Function to compute the Embedded Fractures Preconditioner
function [A] = Compute(obj,A,symMat,varargin)

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   uBlk = obj.getVarBlocks("displacements");
   jBlk = obj.getVarBlocks("fractureJump");

   % If RACP is requested, compute physical stiffness ratio on unscaled blocks
   k_phys = [];
   if strcmpi(obj.regularizationMethod, 'racp')
      k_phys = obj.computePhysicalK(A{uBlk,uBlk}, A{uBlk,jBlk});
   end

   % Compute and apply Ruiz diagonal scaling if requested (ON THE TRUE MATRIX A)
   A = obj.Ruiz.Compute(A);

   if strcmpi(obj.regularizationMethod, 'racp')
      % Physics-weighted RACP augmentation on Ruiz-scaled block
      A22 = obj.augmentA22_RACP(A{uBlk,uBlk}, A{uBlk,jBlk}, A{jBlk,uBlk}, A{jBlk,jBlk}, k_phys);
   else
      % Baseline diagShift: Make A22 slightly more diagonally dominant
      A22 = obj.regularizeA22(A{jBlk,jBlk});
   end

   % Compute the test space
   TV0 = obj.AMG.buildTestSpace(obj.problemsolver, 1:numel(obj.AMG.blockInfo));

   % Scale the test space if needed
   TV0 = obj.Ruiz.applyDinv(TV0, uBlk);

   % Compute factorization of the jump block
   obj.FSAI.Compute(A22, symMat(jBlk,jBlk), true);

   % Compute reverse Schur complement of mechanics block
   invC = @(x) obj.FSAI.ApplyLeft(x);
   obj.S = A{uBlk,uBlk} - A{uBlk,jBlk}*invC(A{jBlk,uBlk});
   fprintf("Condest 11 %e, S %e. Cond 22 %e, 21 %e, 12 %e, 22mod %e\n",...
      condest(A{uBlk,uBlk}),condest(obj.S),cond(full(A{jBlk,jBlk})),cond(full(A{jBlk,uBlk})),cond(full(A{uBlk,jBlk})),cond(full(A22)));

   % Check in case of symmetric indefinite systems for FSAI
   AMGSym = obj.FSAI.PrecSym && obj.PrecSym && obj.FSAI.posDef;

   % Compute the amg for mechanics block
   obj.AMG.Compute(obj.S, AMGSym, TV0, true);

   % Define handles for application of the preconditioner
   obj.Apply_L = @(x) obj.ApplyLeft(x, A{uBlk,jBlk}, A{jBlk,uBlk}, invC);
   obj.Apply_R = @(x) obj.ApplyRight(x);

end


