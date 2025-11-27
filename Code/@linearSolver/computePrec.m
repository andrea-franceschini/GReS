
% Function for the computation of the preconditioner
function computePrec(obj,A)

   if (norm(A-A',"fro")/norm(A,"fro") > 1e-7)
      obj.params.symm = false;
   else
      obj.params.symm = true;
   end

   time_start = tic;
   switch obj.PrecType

      % Compute the AMG preconditioner
      case 'amg'
         
         % Treat Boundary conditions 
         lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
         A(A==1) = lmax/10;

         % Compute the test space
         if(obj.phys == 0) % fluids
            TV0 = ones(size(A,1),1);
         elseif(obj.phys == 1)
            TV0 = mk_rbm_3d(obj.domain.grid.topology.coordinates);
         end
         
         set_DEBINFO();

         % Actually compute the AMG
         obj.Prec = cpt_aspAMG(obj.params,A,TV0);

         % Define Mfun
         obj.MfunL = @(r) AMG_Vcycle(obj.Prec,A,r);
         obj.MfunR = @(r) r;

      % Compute the FSAI preconditioner
      case 'fsai'
         smootherOp = smoother(A,obj.params.symm,obj.params.smoother);

         % Define Mfun
         [obj.MfunL,obj.MfunR] = defineMfunFSAI(obj,smootherOp);
   end
   T_setup = toc(time_start);
   obj.aTimeComp = obj.aTimeComp + T_setup;
   obj.nComp = obj.nComp + 1;
end








% Helper function for computePrec
function set_DEBINFO()
   global DEBINFO;
   % PARTE GENERALE
   DEBINFO.flag = true;
   DEBINFO.flag = false;
   % PARTE PER LA PROLONGATION
   DEBINFO.prol = [];
   % STAMPARE SI/NO
   DEBINFO.prol.prt_flag = false;
   % UNITA DI STAMPA
   DEBINFO.prol.ofile = 0;
   % STAMPA NUMERO DI ITERAZIONI NEL CALCOLO PROL
   DEBINFO.prol.it_print = false;
   % STAMPA LISTA VICINI NEL CALCOLO PROL
   DEBINFO.prol.neigh_print = false;

   % PARTE PER IL COARSENING
   DEBINFO.coarsen = [];
   % STAMPARE SI/NO
   DEBINFO.coarsen.draw_dist = false;
end







% Helper function for computePrec
% Function to determine how MfunL and MfunR are for each fsai preconditioner
function [MfunL,MfunR] = defineMfunFSAI(obj,smootherOp)
   omega = smootherOp.omega;
   if strcmp(lower(obj.params.smoother.method),'afsai_enh')
      F     = smootherOp.left;
      FT    = smootherOp.right;
      W     = smootherOp.W;
      THETA = smootherOp.THETA;
      MfunL = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
      MfunR = @(x) omega*(FT*(F*x) + W*(THETA*(W'*x)));
   else
      if smootherOp.LS_deg > 0
         MfunL = smootherOp.polyPrec;
         MfunR = smootherOp.polyPrec;
      else
         if strcmp(lower(obj.params.smoother.method),'blk_j')
            MfunL = smootherOp.BLKJ;
            MfunR = smootherOp.BLKJ;
         elseif strcmp(lower(obj.params.smoother.method),'bafsai')
            MfunL = smootherOp.BAFSAI;
            MfunR = smootherOp.BAFSAI;
         elseif strcmp(lower(obj.params.smoother.method),'ddsw')
            MfunL = smootherOp.DDSW1;
            MfunR = smootherOp.DDSW2;
         else
            if (numel(smootherOp.left_out) + numel(smootherOp.right_out)) == 0
               % Simple smoother
               MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
               MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
            else
               % Nested smoother
               MfunL = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                   (smootherOp.left_out*(smootherOp.left*x))));
               MfunR = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.right_out*...
                                   (smootherOp.left_out*(smootherOp.left*x))));
            end
         end
      end
   end
end
