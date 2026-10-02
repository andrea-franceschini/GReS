function Compute(obj,A,symm,varargin)

   % Understand if part of a block preconditioner
   if nargin < 4
      block = false;
   else
      block = varargin{1};
   end

   if iscell(A)
      A = A{1,1};
   end

   % If sym == 0 then the matrix is nonsymmetric
   if ~symm
      obj.PrecSym = false;
   else
      obj.PrecSym = true;
   end

   % Treat Boundary conditions if not coming from a block preconditioner 
   if ~block
      warning('off', 'MATLAB:eigs:NotAllEigsConvKeep');
      lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);

      d = diag(A);
      idx = (d == 1);
      d(idx) = lmax/10;
      A = spdiags(d, 0, A);
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