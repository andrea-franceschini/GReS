% Function to compute the Embedded Fractures Preconditioner
function [A] = Compute(obj,A,symMat,varargin)

   if ~iscell(A)
      error("Parent preconditioner not configured to run with multiphys multidomain");
   end

   % Set the symmetry of the full preconditioner
   obj.PrecSym = floor(sum(symMat,'all')/(size(symMat,1))^2);

   uBlk = obj.getVarBlocks("displacements");
   jBlk = obj.getVarBlocks("fractureJump");

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

   % Compute factorization of the jump block
   obj.FSAI.Compute(A{jBlk,jBlk}, symMat(jBlk,jBlk), true);
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


