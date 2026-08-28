function Compute(obj,A,symm,varargin)
   
   % Check inputs are correct 
   if nargin < 4
      gresLog().log(3,'test space not passed to aAMG, using defaults');
      if obj.phys == 0
         TV0 = ones(size(A,1),1);
      elseif obj.phys == 1 
         TV0 = mk_rbm_3d(obj.generalsolver.domains(1).grid.coordinates);
      else
         error('Physics not recognized');
      end
   else
      % Get the test space
      TV0 = varargin{1};
   end

   % Understand if part of a block preconditioner
   if nargin < 5
      block = false;
   else
      block = varargin{2};
   end

   if iscell(A)
      A = A{1,1};
   end

   % If sym == 0 then the matrix is nonsymmetric
   if ~symm
      obj.params.symm = false;
      obj.PrecSym = false;
   else
      obj.params.symm = true;
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

   set_DEBINFO();

   % Compute the AMG preconditioner
   obj.Prec = cpt_aspAMG(obj.params,A,TV0,obj.DEBUGflag);

   % Define Mfun
   obj.Apply_L = @(r) obj.ApplyLeft(r,A);
   obj.Apply_R = @(r) obj.ApplyRight(r);

end








% Helper function for computePrec
function set_DEBINFO()
   global DEBINFO;
   % GENERAL 
   DEBINFO.flag = false;

   % PROLONGATION
   DEBINFO.prol = [];
   % prints
   DEBINFO.prol.prt_flag = false;
   DEBINFO.prol.ofile = 0;
   % iterations in prolongation
   DEBINFO.prol.it_print = false;
   % nearest neighbours print in prol
   DEBINFO.prol.neigh_print = false;

   % COARSENING
   DEBINFO.coarsen = [];
   % prints
   DEBINFO.coarsen.draw_dist = false;
end