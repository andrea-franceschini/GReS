function varargout = printStats(obj, varargin)
% printStats - Display linear solver execution statistics and optionally return diagnostic struct.
%
% Syntax:
%   obj.printStats()
%   obj.printStats('short')
%   obj.printStats('silent')
%   stats = obj.printStats(...)
%
% Options:
%   'short'  - Skip per-solve row table and only print summary metrics
%   'silent' - Suppress all console printing (useful when only extracting stats struct)
%   'ew'     - Display Eisenstat-Walker adaptive tolerance statistics

   % Allocate zero SAM timings if Semi-Algebraic Multigrid was not active
   if ~obj.useSAM
      obj.SAM.CompLin = zeros(size(obj.precCompLin));
   end

   % Parse display flags
   isSilent  = any(strcmpi(varargin, 'silent')) || any(strcmpi(varargin, 'quiet'));
   showShort = any(strcmpi(varargin, 'short'));
   showEW    = any(strcmpi(varargin, 'ew'));

   % Print detailed per-solve diagnostics table unless short or silent display is requested
   if ~showShort && ~isSilent
      fprintf('The preconditioner was computed at time(s):\n');
      for i = 1:numel(obj.whenComputed)
         fprintf('             %d\t%e\n', i, obj.whenComputed(i));
      end

      isRuiz = ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag;

      % Configure table columns based on whether Ruiz scaling was active
      if isRuiz
         sepLine = [repmat('-', 1, 83), '\n'];
         header  = '| %8s | %6s | %4s | %8s | %7s | %8s | %8s | %9s |\n';
         rowFmt  = '| %.2e | %6d | %4d | %.2e | %.1e | %.2e | %8s | %9s |\n';
         colNames = {'PhysTime', 'Sol N.', 'Iter', 'SolTime', 'Symm', 'PrecTime', 'TolReq', 'TolScaled'};
      else
         sepLine = [repmat('-', 1, 71), '\n'];
         header  = '| %8s | %6s | %4s | %8s | %7s | %8s | %8s |\n';
         rowFmt  = '| %.2e | %6d | %4d | %.2e | %.1e | %.2e | %8s |\n';
         colNames = {'PhysTime', 'Sol N.', 'Iter', 'SolTime', 'Symm', 'PrecTime', 'TolReq'};
      end

      fprintf(['\n', sepLine, sprintf(header, colNames{:}), sepLine]);

      % Format each linear solve step row
      for i = 1:size(obj.solveTLin, 2)
         tolReq = NaN;
         if i <= numel(obj.tolReqLin), tolReq = obj.tolReqLin(i); end
         if ~isnan(tolReq), tolReqStr = sprintf('%8.2e', tolReq); else, tolReqStr = '     ---'; end

         rowArgs = {obj.timeLin(i), i, obj.iterLin(i), obj.solveTLin(i), ...
            obj.symFlagLin(i), obj.precCompLin(i) + obj.SAM.CompLin(i), tolReqStr};

         if isRuiz
            tolScaled = NaN;
            if i <= numel(obj.scaledTolLin), tolScaled = obj.scaledTolLin(i); end
            if ~isnan(tolScaled), tolScaledStr = sprintf('%9.2e', tolScaled); else, tolScaledStr = '      ---'; end
            rowArgs{end+1} = tolScaledStr;
         end

         fprintf(rowFmt, rowArgs{:});
      end
      fprintf(sepLine);
   end

   % Print Eisenstat-Walker statistics if requested
   if showEW && ~isSilent && ~isempty(obj.convStrat)
      obj.convStrat.printStats();
   end

   % Print summary statistics
   if ~isSilent
      fprintf('Average Preconditioner computation time = %e\n', obj.aTimeComp / max(1, obj.nComp));
      fprintf('Preconditioner computed %d times\n', numel(obj.whenComputed));
      if obj.useSAM
         nnzSAMcomp = nnz(obj.SAM.CompLin);
         fprintf('Average SAM computation time = %e\n', sum(obj.SAM.CompLin) / max(1, nnzSAMcomp));
         fprintf('SAM computed %d times\n', nnzSAMcomp);
      end
      fprintf('Average Solve time = %e\n', obj.aTimeSolve / max(1, obj.nSolve));
      fprintf('Average number of iterations = %f\n', obj.aIter / max(1, obj.nSolve));
      fprintf('Max number of iterations = %d\n', obj.maxIter);
      fprintf('Total time for computation of the linear systems = %e\n', ...
         obj.aTimeComp + obj.aTimeSolve + sum(obj.SAM.CompLin));
      fprintf('Max system size = %d\n', obj.systemSize);
   end

   % Return structured metrics when requested by caller
   if nargout > 0
      stats = struct();
      stats.nSolve       = obj.nSolve;
      stats.nComp        = numel(obj.whenComputed);
      stats.whenComputed = obj.whenComputed;
      stats.totalIter    = obj.aIter;
      stats.avgIter      = obj.aIter / max(1, obj.nSolve);
      stats.maxIter      = obj.maxIter;
      stats.iterLin      = obj.iterLin;
      stats.timeLin      = obj.timeLin;
      stats.solveTLin    = obj.solveTLin;
      stats.precCompLin  = obj.precCompLin;
      stats.totalTime    = obj.aTimeComp + obj.aTimeSolve + sum(obj.SAM.CompLin);
      stats.systemSize   = obj.systemSize;
      stats.solverType   = obj.SolverType;
      varargout{1} = stats;
   end

end
