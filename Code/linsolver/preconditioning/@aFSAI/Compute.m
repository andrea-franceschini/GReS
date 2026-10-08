function Compute(obj,A,symm,varargin)

   if iscell(A)
      % Ruiz is done only here as it is entered only if the problem to be
      % solved is a single physics using amg as a preconditioner. In the
      % other cases is handled in the block preconditioners

      % Compute and apply Ruiz diagonal scaling if requested
      A = obj.Ruiz.Compute(A);

      % Convert the matrix to sparse double
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
   obj.posDef = obj.Prec.is_posdef;

   % Define Mfun
   obj.Apply_L = @(r) obj.ApplyLeft(r,A);
   obj.Apply_R = @(r) obj.ApplyRight(r);

end