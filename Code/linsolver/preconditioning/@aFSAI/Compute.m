function A = Compute(obj,A,symm,varargin)

   if iscell(A)
      A = A{1,1};
   end

   % If sym == 0 then the matrix is nonsymmetric
   if ~symm
      obj.PrecSym = false;
   else
      obj.PrecSym = true;
   end

   % Compute the FSAI preconditioner
   obj.Prec = smoother(A, symm, obj.param, obj.DEBUGflag);

   % Check if there was the fallback and propagate the symmetry
   if ~obj.Prec.is_posdef && obj.PrecSym == true
      obj.PrecSym = false;
   end

   % Define Mfun
   obj.Apply_L = @(r) obj.ApplyLeft(r,A);
   obj.Apply_R = @(r) obj.ApplyRight(r);

end
