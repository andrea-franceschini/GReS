function varargout = printStats(obj, varargin)
% printStats - Display linear solver execution statistics and optionally return diagnostic struct.
%
% Syntax:
%   obj.printStats()
%   obj.printStats(level)
%   obj.printStats('silent')
%   obj.printStats('short')
%   obj.printStats('default')
%   obj.printStats('detailed')
%   obj.printStats('all')
%   stats = obj.printStats(...)
%
% Verbosity Levels (0 to 4):
%   0 - 'silent', 'quiet'   : Suppress all output.
%   1 - 'short', 'summary'  : Print aggregated summary metrics only.
%   2 - 'default', 'standard': Standard overview table + summary metrics.
%   3 - 'detailed'          : Detailed table (+ NL-It, PrecEff, TolScaled, DeltaT) + summary.
%   4 - 'all', 'long'       : Comprehensive table (+ SAMTime, Symm), EW statistics table (if active).
%
% Additional Options:
%   'ew' - Explicitly display Eisenstat-Walker statistics table regardless of level.

   % Allocate zero SAM timings if Semi-Algebraic Multigrid was not active
   if ~obj.useSAM
      obj.SAM.CompLin = zeros(size(obj.precCompLin));
   end

   % Determine verbosity level (default = 2: 'default' / 'standard')
   verbosity = parseVerbosity(varargin{:});
   showEW = any(strcmpi(varargin, 'ew')) || (verbosity >= 4);

   % Print table if verbosity >= 2
   if verbosity >= 2
      fprintf('The preconditioner was computed at time(s):\n');
      for i = 1:numel(obj.whenComputed)
         fprintf('             %d\t%e\n', i, obj.whenComputed(i));
      end

      % Build modular column definitions
      cols = getTableColumns(obj, verbosity);

      % Build header and separator line
      headerFmt = '|';
      rowFmt    = '|';
      colNames  = cell(1, numel(cols));
      totalWidth = 1; % leading '|'

      for c = 1:numel(cols)
         w = cols(c).width;
         headerFmt = [headerFmt, sprintf(' %%%ds |', w)]; %#ok<AGROW>
         rowFmt    = [rowFmt, sprintf(' %%%ds |', w)];    %#ok<AGROW>
         colNames{c} = cols(c).name;
         totalWidth = totalWidth + w + 3;
      end
      headerFmt = [headerFmt, '\n'];
      rowFmt    = [rowFmt, '\n'];
      sepLine   = [repmat('-', 1, totalWidth), '\n'];

      % Print table header
      fprintf(['\n', sepLine, sprintf(headerFmt, colNames{:}), sepLine]);

      % Format and print each linear solve row
      nSolves = numel(obj.solveTLin);
      for i = 1:nSolves
         rowValues = cell(1, numel(cols));
         for c = 1:numel(cols)
            rowValues{c} = cols(c).getValue(obj, i);
         end
         fprintf(rowFmt, rowValues{:});
      end
      fprintf(sepLine);
   end

   % Print Eisenstat-Walker statistics if active and requested
   if showEW && verbosity > 0 && ~isempty(obj.convStrat)
      obj.convStrat.printStats();
   end

   % Print summary statistics for verbosity >= 1
   if verbosity >= 1
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
      stats.precEffLin   = obj.precEffLin;
      stats.tolReqLin    = obj.tolReqLin;
      stats.scaledTolLin = obj.scaledTolLin;
      stats.tolReachLin  = obj.tolReachLin;
      stats.newtonLin    = obj.newtonLin;
      stats.totalTime    = obj.aTimeComp + obj.aTimeSolve + sum(obj.SAM.CompLin);
      stats.systemSize   = obj.systemSize;
      stats.solverType   = obj.SolverType;
      varargout{1} = stats;
   end

end

% =========================================================================
% Helper: Parse verbosity level from arguments
% =========================================================================
function verbosity = parseVerbosity(varargin)
   verbosity = 2; % default ('default' / 'standard')

   for i = 1:numel(varargin)
      arg = varargin{i};
      if isnumeric(arg) && isscalar(arg)
         verbosity = max(0, min(4, round(arg)));
         return;
      elseif ischar(arg) || isstring(arg)
         switch lower(char(arg))
            case {'silent', 'quiet'}
               verbosity = 0;
               return;
            case {'short', 'summary'}
               verbosity = 1;
               return;
            case {'default', 'standard'}
               verbosity = 2;
            case {'detailed'}
               verbosity = 3;
            case {'all', 'long'}
               verbosity = 4;
         end
      end
   end
end

% =========================================================================
% Helper: Build active modular table columns based on verbosity and availability
% =========================================================================
function activeCols = getTableColumns(obj, verbosity)

   isRuiz = ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag;
   hasRecycledPrec = ~isempty(obj.precEffLin) && any(obj.precEffLin > 0 & obj.precEffLin ~= obj.precCompLin);
   hasNewton = ~isempty(obj.newtonLin) && any(obj.newtonLin > 0);
   hasDeltaT = ~isempty(obj.Delta_T) && any(obj.Delta_T > 0);

   % Column Registry:
   % name: header display title
   % width: field width in formatted table
   % minVerb: minimum verbosity level (2: default, 3: detailed, 4: all)
   % isAvail: predicate checking whether column has data
   % getValue: closure returning formatted string for row i
   allCols = [ ...
      struct('name', 'PhysTime', 'width', 8, 'minVerb', 2, 'isAvail', @() true, ...
             'getValue', @(o, i) formatNum('%.2e', o.timeLin, i, 8)), ...
      struct('name', 'Sol N.', 'width', 6, 'minVerb', 2, 'isAvail', @() true, ...
             'getValue', @(o, i) formatInt('%6d', i)), ...
      struct('name', 'NL-It', 'width', 5, 'minVerb', 3, 'isAvail', @() hasNewton, ...
             'getValue', @(o, i) formatNum('%5d', o.newtonLin, i, 5)), ...
      struct('name', 'Iter', 'width', 4, 'minVerb', 2, 'isAvail', @() true, ...
             'getValue', @(o, i) formatNum('%4d', o.iterLin, i, 4)), ...
      struct('name', 'SolTime', 'width', 8, 'minVerb', 2, 'isAvail', @() true, ...
             'getValue', @(o, i) formatNum('%.2e', o.solveTLin, i, 8)), ...
      struct('name', 'PrecTime', 'width', 8, 'minVerb', 2, 'isAvail', @() true, ...
             'getValue', @(o, i) formatPrecTime(o, i)), ...
      struct('name', 'PrecEff', 'width', 8, 'minVerb', 3, 'isAvail', @() hasRecycledPrec, ...
             'getValue', @(o, i) formatNum('%.2e', o.precEffLin, i, 8)), ...
      struct('name', 'TolReq', 'width', 8, 'minVerb', 2, 'isAvail', @() true, ...
             'getValue', @(o, i) formatNum('%.2e', o.tolReqLin, i, 8)), ...
      struct('name', 'TolReach', 'width', 8, 'minVerb', 2, 'isAvail', @() true, ...
             'getValue', @(o, i) formatNum('%.2e', o.tolReachLin, i, 8)), ...
      struct('name', 'TolScaled', 'width', 9, 'minVerb', 3, 'isAvail', @() isRuiz, ...
             'getValue', @(o, i) formatNum('%.2e', o.scaledTolLin, i, 9)), ...
      struct('name', 'DeltaT', 'width', 8, 'minVerb', 3, 'isAvail', @() hasDeltaT, ...
             'getValue', @(o, i) formatNum('%.2e', o.Delta_T, i, 8)), ...
      struct('name', 'SAMTime', 'width', 8, 'minVerb', 4, 'isAvail', @() obj.useSAM, ...
             'getValue', @(o, i) formatNum('%.2e', o.SAM.CompLin, i, 8)), ...
      struct('name', 'Symm', 'width', 7, 'minVerb', 4, 'isAvail', @() true, ...
             'getValue', @(o, i) formatNum('%.1e', o.symFlagLin, i, 7)) ...
   ];

   % Filter active columns based on verbosity and availability
   keep = false(1, numel(allCols));
   for c = 1:numel(allCols)
      if verbosity >= allCols(c).minVerb && allCols(c).isAvail()
         keep(c) = true;
      end
   end
   activeCols = allCols(keep);
end

% =========================================================================
% Formatting helper closures
% =========================================================================
function str = formatNum(fmt, arr, idx, width)
   if idx <= numel(arr) && ~isnan(arr(idx))
      str = sprintf(fmt, arr(idx));
      if length(str) < width
         str = [repmat(' ', 1, width - length(str)), str];
      end
   else
      str = [repmat(' ', 1, max(0, width - 3)), '---'];
   end
end

function str = formatInt(fmt, val)
   str = sprintf(fmt, val);
end

function str = formatPrecTime(obj, idx)
   val = 0;
   if idx <= numel(obj.precCompLin)
      val = obj.precCompLin(idx);
   end
   if obj.useSAM && idx <= numel(obj.SAM.CompLin)
      val = val + obj.SAM.CompLin(idx);
   end
   str = sprintf('%.2e', val);
   if length(str) < 8
      str = [repmat(' ', 1, 8 - length(str)), str];
   end
end
