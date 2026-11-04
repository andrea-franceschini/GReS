function printStats(obj, varargin)
   % Allocate for stats if sam not used
   if ~obj.useSAM
      obj.SAM.CompLin = zeros(size(obj.precCompLin));
   end
   
   % Parse optional arguments
   showShort = any(strcmpi(varargin, 'short'));
   showEW    = any(strcmpi(varargin, 'ew'));
   
   % If not 'short', show full linear solver table and preconditioner times
   if ~showShort
      fprintf('The preconditioner was computed at time(s):\n');
      for i = 1:numel(obj.whenComputed)
         fprintf('             %d\t%e\n', i, obj.whenComputed(i));
      end
   
      isRuiz = ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag;
   
      % Setup table headers and formats dynamically
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
   
      % Helper to format optional tolerance values or fallback dashes
      fmtTol = @(val, width, fmt) ternary(~isnan(val), sprintf(fmt, val), sprintf(['%', width, 's'], '---'));
   
      for i = 1:size(obj.solveTLin, 2)
         % Format TolReq
         tolReq = NaN;
         if i <= numel(obj.tolReqLin), tolReq = obj.tolReqLin(i); end
         tolReqStr = ternary(~isnan(tolReq), sprintf('%8.2e', tolReq), '     ---');
   
         rowArgs = {obj.timeLin(i), i, obj.iterLin(i), obj.solveTLin(i), ...
            obj.symFlagLin(i), obj.precCompLin(i) + obj.SAM.CompLin(i), tolReqStr};
   
         if isRuiz
            tolScaled = NaN;
            if i <= numel(obj.scaledTolLin), tolScaled = obj.scaledTolLin(i); end
            rowArgs{end+1} = ternary(~isnan(tolScaled), sprintf('%9.2e', tolScaled), '      ---');
         end
   
         fprintf(rowFmt, rowArgs{:});
      end
      fprintf(sepLine);
   end
   
   % Print Eisenstat-Walker statistics if explicitly requested
   if showEW && ~isempty(obj.convStrat)
      obj.convStrat.printStats();
   end
   
   % Summary statistics
   fprintf('Average Preconditioner computation time = %e\n', obj.aTimeComp / obj.nComp);
   fprintf('Preconditioner computed %d times\n', numel(obj.whenComputed));
   if obj.useSAM
      nnzSAMcomp = nnz(obj.SAM.CompLin);
      fprintf('Average SAM computation time = %e\n', sum(obj.SAM.CompLin) / nnzSAMcomp);
      fprintf('SAM computed %d times\n', nnzSAMcomp);
   end
   fprintf('Average Solve time = %e\n', obj.aTimeSolve / obj.nSolve);
   fprintf('Average number of iterations = %f\n', obj.aIter / obj.nSolve);
   fprintf('Max number of iterations = %d\n', obj.maxIter);
   fprintf('Total time for computation of the linear systems = %e\n', ...
      obj.aTimeComp + obj.aTimeSolve + sum(obj.SAM.CompLin));
   fprintf('Max system size = %d\n', obj.systemSize);
   end
   
   % Local inline ternary helper
   function out = ternary(cond, a, b)
   if cond, out = a; else, out = b; end
end