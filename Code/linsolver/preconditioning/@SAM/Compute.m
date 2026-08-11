function Compute(obj,A,symm,varargin)
   % Get input preconditioner
   Prec = varargin{1};

   if obj.useSAM
      if isempty(obj.N) && ~obj.requestComp
         % If it is empty it was not computed then set to zero the setup
         % time
         obj.tSetup = 0;
         
         % Set up the usage when the sam is the identity
         obj.Apply_L = @(x) Prec.Apply_L(x);
         obj.Apply_R = @(x) Prec.Apply_R(x);

      elseif obj.requestComp
         % Compute the SAM
         time_start = tic;
         [obj.N, ~] = MEX_sam_adaptive_left(A,obj.A0,obj.maxThreads,obj.nstep,obj.stp_size,obj.epss);
      
         % Symmetrize
         if obj.PrecSym && symm
            obj.N = 0.5*(obj.N+obj.N');
         end
         obj.tSetup = toc(time_start);
   
         % Take the stats
         obj.aTimeComp = obj.aTimeComp + obj.tSetup;
         obj.nComp = obj.nComp + 1;
         obj.nSolveSinceLastComp = 0;

         % Set up the usage 
         obj.Apply_L = @(x) obj.ApplyLeft( x,Prec);
         obj.Apply_R = @(x) obj.ApplyRight(x,Prec);
   
      else
         % Reusing the SAM, just add one to solves since the computation
         obj.nSolveSinceLastComp = obj.nSolveSinceLastComp + 1;
         obj.tSetup = 0;
      end
   else
      obj.tSetup = 0;

      % Set up the usage when the sam is not requested
      obj.Apply_L = @(x) Prec.Apply_L(x);
      obj.Apply_R = @(x) Prec.Apply_R(x);
   end
end