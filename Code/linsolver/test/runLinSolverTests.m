function allPassed = runLinSolverTests(varargin)
% runLinSolverTests - Test suite and regression verification for linearSolver.
%
% Usage:
%   runLinSolverTests()
%   runLinSolverTests('includeSimulation', true)
%   runLinSolverTests('updateReference', true)

p = inputParser;
addParameter(p, 'updateReference', false, @islogical);
addParameter(p, 'includeSimulation', true, @islogical);
parse(p, varargin{:});

updateReference   = p.Results.updateReference;
includeSimulation = p.Results.includeSimulation;

origDir = pwd;
testDir = fileparts(mfilename('fullpath'));
cd(testDir);
cleanUpObj = onCleanup(@() cd(origDir));

% Set for reproducibility
rng(0);

refFile = fullfile(testDir, 'matrices', 'linsolver_test_reference.mat');

% Test case definitions covering single-phase flow, mechanics, coupled Biot,
% stick-slip contact, embedded fracture, two-phase flow, and multiphase poromechanics.
testCases = {
   struct('file', 'matrices/spflow_fluxBarrier.mat', 'rhsFlag', 0, 'desc', 'Single-Phase Flow (flux barrier)');
   struct('file', 'matrices/mech_treporti.mat',      'rhsFlag', 0, 'desc', 'Linear Mechanics (Treporti)');
   struct('file', 'matrices/richards_column.mat',    'rhsFlag', 0, 'desc', 'Richards Equation (column)');
   struct('file', 'matrices/biotFEM_treporti.mat',   'rhsFlag', 0, 'desc', 'Coupled Biot FEM (surface load)');
   struct('file', 'matrices/stickslipEFEM.mat',      'rhsFlag', 0, 'desc', 'Embedded Fracture (EFEM)');
   struct('file', 'matrices/twoPhase_SPE.mat',        'rhsFlag', 0, 'desc', 'Two-Phase Flow (SPE)');
   struct('file', 'matrices/mpPoro_staircase.mat',    'rhsFlag', 0, 'desc', 'Multiphase Poromechanics (3-field)');
   struct('file', 'matrices/spflow_fluxBarrier.mat', 'rhsFlag', 1, 'desc', 'Flow with b = ones(n, 1)');
   struct('file', 'matrices/spflow_fluxBarrier.mat', 'rhsFlag', 2, 'desc', 'Flow with b = rand(n, 1) [seeded]');
   struct('file', 'matrices/mech_treporti.mat',      'rhsFlag', 3, 'desc', 'Mechanics with manufactured A*ones');
   struct('file', 'matrices/biotFEM_treporti.mat',   'rhsFlag', 3, 'desc', 'Biot with manufactured A*ones');
   struct('file', 'matrices/twoPhase_SPE.mat',        'rhsFlag', 1, 'desc', 'Two-Phase Flow with b = ones')
   };

nTests = numel(testCases);

% Load reference baseline if existing, otherwise force baseline creation
refData = [];
if ~updateReference && isfile(refFile)
   loaded = load(refFile);
   if isfield(loaded, 'referenceResults')
      refData = loaded.referenceResults;
   end
else
   updateReference = true;
end

results = cell(nTests, 1);
reports = cell(nTests, 1);
allPassed = true;

% Print run mode header
if updateReference
   fprintf('Updating reference baseline\n\n');
else
   fprintf('Running regression verification\n\n');
end

% Execute each standalone matrix test case sequentially
for k = 1:nTests
   tc = testCases{k};
   fprintf('[%2d/%2d] Running: %-32s (rhsFlag=%d)... ', k, nTests, tc.file, tc.rhsFlag);

   % Seed random number generator to ensure reproducible RHS vectors for random test cases
   if tc.rhsFlag == 2
      rng(0);
   end

   % Solve test case via solveSingle
   [x, stats, flag, relres, iter, resvec] = linearSolver.solveSingle(tc.file, 'rhsFlag', tc.rhsFlag, 'verb', false);

   % Fail fast on divergence
   assert(flag == 0, 'Test %s (rhsFlag=%d) failed to converge!', tc.file, tc.rhsFlag);
   fprintf('Converged (%d iters, relres=%.2e)\n', iter, relres);

   % Package full solution result struct
   res = stats;
   res.file = tc.file;
   res.rhsFlag = tc.rhsFlag;
   res.desc = tc.desc;
   res.x = x;
   results{k} = res;

   % Evaluate comparison against reference baseline
   refItem = [];
   if ~updateReference && k <= numel(refData)
      refItem = refData{k};
   end
   [rep, passed] = evaluateReport(tc, x, stats, relres, iter, refItem);
   reports{k} = rep;
   allPassed = allPassed && passed;
end

% Execute nonlinear simulations if enabled
simReport = [];
sedimReport = [];
if includeSimulation
   fprintf('\n[SIM ] Running: Stick-Slip contact simulation (simul/runStickSlipOpen.m)...\n');
   run('simul/runStickSlipOpen.m');

   % Extract simulation linear solver performance metrics via printStats
   simStats = solver.linsolver.printStats('short');
   simPassed = (solver.t >= 8.0) && isa(solver.linsolver.Prec, 'RACP') && isa(solver.linsolver.convStrat, 'EisenstatWalker');
   allPassed = allPassed && simPassed;
   simReport = createSimulationReport(simStats, simPassed, 'RACP');

   % Execute evolving grid sedimentation simulation to verify the growing preconditioner
   fprintf('\n[SIM ] Running: Sedimentation 2D simulation (simul/runLinearSedim.m)...\n');
   run('simul/runLinearSedim.m');

   % Verify only that the sedimentation simulation completes successfully with growing preconditioner
   sedimStats = solverSedim.linsolver.printStats('short');
   sedimPassed = (solverSedim.t >= 40.0) && isa(solverSedim.linsolver.Prec, 'growing');
   allPassed = allPassed && sedimPassed;
   sedimReport = createSimulationReport(sedimStats, sedimPassed, 'growing');
end

% Persist updated reference baseline to MAT file
if updateReference
   referenceResults = results;
   if includeSimulation
      referenceSimStats = simStats;
      save(refFile, 'referenceResults', 'referenceSimStats');
   else
      save(refFile, 'referenceResults');
   end
   fprintf('\nReference baseline successfully saved to: %s\n', refFile);
end

% Render summary table and enforce regression assertions
printSummaryTable(reports, simReport, sedimReport, updateReference);
assert(allPassed, 'One or more tests regressed or drifted from reference baseline!');

end


% =========================================================================
% Local Helper Functions
% =========================================================================

% Evaluate test metrics and compare against baseline reference if available
function [rep, passed] = evaluateReport(tc, x, stats, relres, iter, refItem)
   rep = struct();
   rep.file = tc.file;
   rep.rhsFlag = tc.rhsFlag;
   rep.solverType = stats.solverType;
   rep.preconditioner = stats.preconditioner;
   rep.iter = iter;
   rep.relres = relres;
   rep.trueResNorm = stats.trueResNorm;
   rep.precTime = stats.precTime;
   rep.solveTime = stats.solveTime;
   
   % Baseline generation mode (no reference comparison needed)
   if isempty(refItem)
      rep.hasRef = false;
      rep.refIter = iter;
      rep.refSolveTime = stats.solveTime;
      rep.relDiffX = 0;
      rep.status = 'SAVED';
      passed = true;
      return;
   end
   
   % Regression verification against saved baseline
   rep.hasRef = true;
   rep.refIter = refItem.iter;
   rep.refSolveTime = refItem.solveTime;
   
   % Measure relative solution vector difference: ||x - xref|| / ||xref||
   rep.relDiffX = norm(x - refItem.x) / max(norm(refItem.x), eps);
   
   % Verify against tolerances. A 10% iteration window accounts for minor
   % platform or BLAS scheduling variations, while 10*tol detects numerical drift.
   maxAllowedDiffX = max(1e-4, 10 * stats.tol);
   iterDelta = iter - refItem.iter;
   iterTol = max(5, ceil(0.1 * refItem.iter));
   
   passed = true;
   rep.status = 'MATCH';
   
   if rep.relDiffX > maxAllowedDiffX
      rep.status = 'DRIFT';
      passed = false;
   elseif iterDelta > iterTol
      rep.status = 'REGRESSION';
      passed = false;
   elseif iterDelta < -iterTol
      rep.status = 'IMPROVED';
   end
end

% Package simulation report structure
function rep = createSimulationReport(stats, passed, precName)
   rep = struct();
   rep.nSolve = stats.nSolve;
   rep.solverType = stats.solverType;
   rep.preconditioner = precName;
   rep.totalIter = stats.totalIter;
   rep.nComp = stats.nComp;
   rep.totalTime = stats.totalTime;
   if passed
      rep.status = 'PASSED';
   else
      rep.status = 'FAILED';
   end
end

% Print clean summary table displaying execution metrics and regression status
function printSummaryTable(reports, simReport, sedimReport, isUpdateMode)
   fprintf('\n');
   if isUpdateMode
      fprintf(' #  %-29s RHS  %-7s  %-13s  Iter    RelRes     TrueResNorm   Status\n', ...
         'Matrix File', 'Solver', 'Precond');
      fprintf('%s\n', repmat('-', 1, 102));
   else
      fprintf(' #  %-29s RHS  %-7s  %-13s  Iter(Ref)   RelRes     TrueResNorm   ||x-xref||/||x||  Status\n', ...
         'Matrix File', 'Solver', 'Precond');
      fprintf('%s\n', repmat('-', 1, 124));
   end
   
   for k = 1:numel(reports)
      r = reports{k};
      [~, baseName, ext] = fileparts(r.file);
      dispName = [baseName, ext];
      if length(dispName) > 29
         dispName = [dispName(1:26), '...'];
      end
   
      precShort = r.preconditioner;
      if length(precShort) > 12
         precShort = [precShort(1:9), '...'];
      end
   
      if isUpdateMode
         fprintf('%2d  %-29s [%d]  %-7s  %-12s %4d      %.2e   %.2e      %-8s\n', ...
            k, dispName, r.rhsFlag, r.solverType, precShort, ...
            r.iter, r.relres, r.trueResNorm, r.status);
      else
         iterStr = sprintf('%3d (%3d)', r.iter, r.refIter);
         fprintf('%2d  %-29s [%d]  %-7s  %-12s  %9s    %.2e   %.2e      %.2e          %-8s\n', ...
            k, dispName, r.rhsFlag, r.solverType, precShort, ...
            iterStr, r.relres, r.trueResNorm, r.relDiffX, r.status);
      end
   end
   
   % Print Stick-Slip contact simulation row if executed
   if ~isempty(simReport)
      if isUpdateMode
         fprintf('%s\n', repmat('-', 1, 102));
         fprintf('SIM %-29s [%d] %-7s  %-12s   %4d    %2d prec  %6.2fs total   %-8s\n', ...
            'StickSlipOpen (Simulation)', simReport.nSolve, upper(simReport.solverType), simReport.preconditioner, ...
            simReport.totalIter, simReport.nComp, simReport.totalTime, simReport.status);
      else
         fprintf('%s\n', repmat('-', 1, 124));
         fprintf('SIM %-29s [%d] %-7s  %-12s   %4d       %2d  prec  %6.2fs total     ----            %-8s\n', ...
            'StickSlipOpen (Simulation)', simReport.nSolve, upper(simReport.solverType), simReport.preconditioner, ...
            simReport.totalIter, simReport.nComp, simReport.totalTime, simReport.status);
      end
   end
   
   % Print Sedimentation growing preconditioner simulation row if executed
   if ~isempty(sedimReport)
      if isUpdateMode
         fprintf('SIM %-29s [%d] %-7s  %-12s   %4d   %2d  prec  %6.2fs total   %-8s\n', ...
            'LinearSedim (Growing Prec)', sedimReport.nSolve, upper(sedimReport.solverType), sedimReport.preconditioner, ...
            sedimReport.totalIter, sedimReport.nComp, sedimReport.totalTime, sedimReport.status);
      else
         fprintf('SIM %-29s [%d] %-7s  %-12s   %4d       %2d  prec  %6.2fs total     ----            %-8s\n', ...
            'LinearSedim (Growing Prec)', sedimReport.nSolve, upper(sedimReport.solverType), sedimReport.preconditioner, ...
            sedimReport.totalIter, sedimReport.nComp, sedimReport.totalTime, sedimReport.status);
      end
   end
   fprintf('\n');
end
