% Function to compute the Embedded Fractures Preconditioner
function [A] = Compute(obj,A,symMat,varargin)

   if ~iscell(A)
      error("Parent preconditioner not configured to run with multiphys multidomain");
   end

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   uBlk = obj.getVarBlocks("displacements");
   jBlk = obj.getVarBlocks("fractureJump");

   % If RACP is requested, compute physical stiffness ratio on unscaled blocks
   k_phys = [];
   if strcmpi(obj.regularizationMethod, 'racp')
      k_phys = obj.computePhysicalK(A{uBlk,uBlk}, A{uBlk,jBlk});
   end

   % Compute or retrieve test space
   if nargin >= 4 && ~isempty(varargin{1})
      TV0 = varargin{1};
      TV0 = TV0(uBlk);
   else
      TV0 = obj.getBlockTestSpace(obj.problemsolver);
      TV0 = TV0(uBlk);
   end

   % Compute and apply Ruiz diagonal scaling if requested
   A = obj.Ruiz.Compute(A);
   TV0 = obj.Ruiz.applyDinv(TV0, uBlk);

   if strcmpi(obj.regularizationMethod, 'racp')
      % Physics-weighted RACP augmentation on Ruiz-scaled block
      A22 = obj.augmentA22_RACP(A{uBlk,uBlk}, A{uBlk,jBlk}, A{jBlk,uBlk}, A{jBlk,jBlk}, k_phys);
   else
      % Baseline diagShift: Make A22 slightly more diagonally dominant
      A22 = obj.regularizeA22(A{jBlk,jBlk});
   end

   % Compute factorization of the jump block
   obj.FSAI.Compute(A22, symMat(jBlk,jBlk), true);
   invC = @(x) obj.FSAI.ApplyLeft(x);

   % Compute reverse Schur complement of mechanics block
   obj.S = A{uBlk,uBlk} - A{uBlk,jBlk}*invC(A{jBlk,uBlk});
   
   % Check in case of symmetric indefinite systems for FSAI
   innerPrecSym = obj.FSAI.PrecSym && obj.PrecSym && obj.FSAI.posDef;

   % Compute the innerPrec for mechanics block
   obj.innerPrec.Compute(obj.S, innerPrecSym, TV0, true);

   % Define handles for application of the preconditioner
   obj.Apply_L = @(x) obj.ApplyLeft(x, A{uBlk,jBlk}, A{jBlk,uBlk}, invC);
   obj.Apply_R = @(x) obj.ApplyRight(x);

end


