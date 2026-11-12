function A = Compute(obj,A,symm,varargin)
   
   % Condense the domains
   A = obj.condenseDomains(A);

   % Check inputs are correct 
   if nargin < 4
      gresLog().log(3,'test space not passed to multiDomAMG, using defaults\n');
      TV0 = obj.buildTestSpace(obj.generalsolver, 1:numel(obj.blockInfo));
   else
      % Get the test space
      TV0 = varargin{1};
   end

   if iscell(A)
      % Ruiz is done only here as it is entered only if the problem to be
      % solved is a single physics using amg as a preconditioner. In the
      % other cases is handled in the block preconditioners

      % Compute and apply Ruiz diagonal scaling if requested
      if ~isempty(obj.Ruiz)
         A = obj.Ruiz.Compute(A);
         TV0 = obj.Ruiz.applyDinv(TV0);
      end

      % Convert the matrix to sparse double
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

   set_DEBINFO();

   % Compute the AMG preconditioner
   obj.Prec = cpt_aspAMG(obj.params,A,TV0,obj.DEBUGflag);
   
   % Find out if it is an indefinite case
   obj.posDef = obj.Prec.isPosDef;

   % Get AMG hierarchy information
   obj.AMG_info = get_AMG_info(obj.Prec,A);

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
