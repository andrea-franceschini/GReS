function computeSinglePhPrec(obj,A,symMat)
   
   if iscell(A)
      A = A{1,1};
   end

   % If symMat == 0 then the matrix is nonsymmetric
   if ~symMat
      obj.params.symm = false;
   else
      obj.params.symm = true;
   end

   switch obj.PrecType

      % Compute the AMG preconditioner
      case 'amg'
         
         % Treat Boundary conditions 
         warning('off', 'MATLAB:eigs:NotAllEigsConvKeep');
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
            if obj.DEBUGflag
               obj.TV0 = TV0;
            end
         end

         % Apply the Ruiz scaling factor to the test space
         if ~isempty(obj.D)
            D = vertcat(obj.D{1:obj.nDom});
            TV0 = diag(D)\TV0;
         end

         set_DEBINFO();

         % Actually compute the AMG
         obj.Prec = cpt_aspAMG(obj.params,A,TV0,obj.DEBUGflag);

         % Define Mfun
         obj.Apply_L = @(r) AMG_Vcycle(obj.Prec,A,r);
         obj.Apply_R = @(r) r;

      % Compute the FSAI preconditioner
      case 'fsai'
         smootherOp = smoother(A,obj.params.symm,obj.params.smoother);

         % Define Mfun
         [obj.Apply_L,obj.Apply_R] = defineMfunFSAI(obj,smootherOp);
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
   % Simple smoother
   MfunL = @(x) omega*(smootherOp.right*(smootherOp.left*x));
   MfunR = @(x) omega*(smootherOp.right*(smootherOp.left*x));
end
