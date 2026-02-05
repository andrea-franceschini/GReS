function computeSinglePhPrec(obj,A)
   
   if iscell(A)
      A = A{1,1};
   end

   if obj.DEBUGflag
      fprintf('\nsymmetry = %e\n\n',norm(A-A','f')/norm(A,'f'));
   end
   if (norm(A-A','f')/norm(A,'f') > obj.nsyTol)
      obj.params.symm = false;
      if obj.DEBUGflag
         fprintf('matrix nonsymmatric %e\n',norm(A-A','f')/norm(A,'f'));
      end
   else
      obj.params.symm = true;
   end

   switch obj.PrecType

      % Compute the AMG preconditioner
      case 'amg'
         
         % Treat Boundary conditions 
         lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);

         d = diag(A);
         idx = (d == 1);
         d(idx) = lmax/10;
         A = spdiags(d, 0, A);

         % Compute the test space
         if(obj.phys == 0) % fluids
            TV0 = ones(size(A,1),1);
         elseif(obj.phys == 1 || obj.phys == 1.1) % true contact mechanichs physics is 1.1, general poromechanics is 1
            TV0 = [];
            for i = 1:obj.nDom
               TV = mk_rbm_3d(obj.domain(i).grid.topology.coordinates);
               TV0 = [TV0;TV];
            end
         end

         
         % coord = obj.domain.grid.topology.coordinates;
         % save("TV0.mat","TV0");
         %save("mat_new.mat","A","TV0","coord");
         % error('ciao');
         set_DEBINFO();

         % Actually compute the AMG
         obj.Prec = cpt_aspAMG(obj.params,A,TV0,obj.DEBUGflag);

         % Define Mfun
         obj.MfunL = @(r) AMG_Vcycle(obj.Prec,A,r);
         obj.MfunR = @(r) r;

      % Compute the FSAI preconditioner
      case 'fsai'
         smootherOp = smoother(A,obj.params.symm,obj.params.smoother);

         % Define Mfun
         [obj.MfunL,obj.MfunR] = defineMfunFSAI(obj,smootherOp);
      otherwise
         error('Non defined preconditioner case')
   end
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
