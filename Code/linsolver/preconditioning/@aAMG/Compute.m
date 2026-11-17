function A = Compute(obj,A,symm,varargin)
   
   % Check inputs are correct 
   if nargin < 4 || isempty(varargin{1})
      TV0 = obj.getBlockTestSpace(obj.generalsolver);
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
      A = cell2matrix(A);
   end

   % Extract assembled numeric test space matrix for cpt_aspAMG
   TV0 = preconditioner.extractTestSpaceMatrix(TV0);

   % If sym == 0 then the matrix is nonsymmetric
   if ~symm
      obj.params.symm = false;
      obj.PrecSym = false;
   else
      obj.params.symm = true;
      obj.PrecSym = true;
   end

   resetTimersAndDebinfo();

   % Compute the AMG preconditioner
   tTotal = tic;
   obj.Prec = cpt_aspAMG(obj.params,A,TV0,obj.DEBUGflag);
   totalTime = toc(tTotal);

   % Record hierarchy construction times for this compute call
   obj.recordBuildTimes(totalTime);

   % Find out if it is an indefinite case
   obj.posDef = obj.Prec.isPosDef;

   % Get AMG hierarchy information
   obj.AMG_info = get_AMG_info(obj.Prec,A);

   % Define Mfun
   obj.Apply_L = @(r) obj.ApplyLeft(r,A);
   obj.Apply_R = @(r) obj.ApplyRight(r);

end








% Helper function to reset debinfo and AMG hierarchy construction timers
function resetTimersAndDebinfo()
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

   % Reset AMG construction timers
   amg_times;
   T_smoo = 0; T_tspa = 0; T_coar = 0; T_CR = 0;
   T_prol = 0; T_FilProl = 0; T_FilCLev = 0; T_LowRank = 0; T_Total = 0;
end
