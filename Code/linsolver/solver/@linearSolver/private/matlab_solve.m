function [x,flag] = matlab_solve(obj,A,b,time,nonlinIter)
% MATLAB_SOLVE Direct backslash solver fallback for linearSolver

   gresLog().log(4,'Fallback to matlab due to size or chronos inexistance\n');

   % Get the full rhs
   if iscell(b)
      b = -cell2matrix(b);
   end

   startT = tic;
   % Solve the system
   A = cell2matrix(A);
   x = A\b;
   Tend = toc(startT);

   % Get the size of the system
   obj.systemSize = max(size(A,1),obj.systemSize);

   symValue = norm(A-A','f')/norm(A,'f');
   obj.params.iter = 0;
   obj.params.actualSolver = 'direct';
   fillStats(obj,Tend,time,symValue,0,0,0,NaN,NaN,nonlinIter,0);
   obj.Delta_T(obj.nSolve) = 0;
   flag = 0;
end
