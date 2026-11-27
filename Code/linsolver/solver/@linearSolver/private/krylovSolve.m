function [x,flag,iter] = krylovSolve(obj,Amat,b)
% KRYLOVSOLVE Dispatches iterative solve to GMRES or SQMR

   switch obj.SolverType
      case 'gmres'
   
         % Solve the system by GMRES
         [x,flag,obj.params.lastRelres,iter1,obj.params.resvec] = gmres_RIGHT(Amat,b,obj.params.restart,obj.toleranceStrategy.Tol,...
            obj.params.maxit/obj.params.restart,...
            obj.SAM.Apply_L,obj.SAM.Apply_R,obj.x0,obj.DEBUGflag);
         iter = (iter1(1) - 1) * obj.params.restart + iter1(2);
   
      case 'sqmr'
   
         % Solve the system by SQMR
         Afun = @(x) Amat*x;
         [x,flag,obj.params.lastRelres,iter,obj.params.resvec] = SQMR(Afun,b,obj.toleranceStrategy.Tol,obj.params.maxit,...
            obj.SAM.Apply_L,obj.SAM.Apply_R,obj.x0,obj.DEBUGflag);
   end

   % Sanity check for all real solutions
   if any(~isreal(x),'all')
      gresLog().warning(3,'wtf')
      x = real(x);
   end
end
