% Save current linear system matrix, RHS, grid coordinates, and solver metadata to .mat file.
% This enables standalone reproduction and benchmark validation of linear systems via solveSingle.
function saveLinearSystem(obj, A, b, currentSolve, varargin)

% Determine if matrix export is requested for this solve step
forcing = ~isempty(varargin) && (isequal(varargin{1}, true) || strcmpi(string(varargin{1}), "true"));
shouldSave = forcing || strcmpi('all', obj.saveMatType) || ...
   (strcmpi('selected', obj.saveMatType) && ismember(currentSolve, obj.saveMatSolves));
if ~shouldSave
   return;
end

% Extract spatial grid coordinates from domain models
coordinates = getCoordinates(obj.generalsolver);
filename = sprintf('linearSystem_solve%d.mat', currentSolve);

% Extract scaling flag and active preconditioner class
ruizFlag = ~isempty(obj.Ruiz) && obj.Ruiz.scalingFlag;
if ~isempty(obj.Prec), preconditioner = class(obj.Prec); else, preconditioner = ''; end

% For fixed-stress split systems, precompute and record relaxation matrix R.
% The relaxation matrix (e.g. volumetric fluid coupling in Biot or multiphase poromechanics)
% is required by fixedStress to decouple mechanics and flow without rebuilding domain matrices.
timeRelaxation = 0;
R = [];
props = struct();
if isa(obj.Prec, 'fixedStress')
   dom1 = obj.generalsolver.domains(1);
   if ismember("MultiphasePoromechanics", dom1.solverNames)
      biot = dom1.getPhysicsSolver("MultiphasePoromechanics");
   elseif ismember("BiotFullyCoupled", dom1.solverNames)
      biot = dom1.getPhysicsSolver("BiotFullyCoupled");
   else
      error('Configuration not recognized for fixed stress in saveLinearSystem');
   end

   if ~isempty(biot)
      t0_R = tic;
      biot.computeRelaxationMatrix();
      timeRelaxation = toc(t0_R);
      R = biot.R;

      % For multiphase poromechanics, extract fluid saturation and phase densities
      % to enable exact Schur complement flow scaling in standalone solveSingle runs
      if ismember("MultiphasePoromechanics", dom1.solverNames)
         pState = dom1.getState("pressure");
         swState = dom1.getState("saturation");
         cIdx = dom1.dofm.getFieldCells("pressure");
         tpSolver = biot.getFlowSolver();
         sp = pState(cIdx);
         props.sw = swState(cIdx);
         props.rhoW = tpSolver.props.wet.density(sp);
         props.rhoNW = tpSolver.props.nwet.density(sp);
      end
   end
end

% Extract physical variable names from domain degree of freedom managers
if ~isempty(obj.generalsolver.domains)
   physCell = arrayfun(@(x) x.dofm.getVariableNames(), obj.generalsolver.domains, 'UniformOutput', false);
   physics = [physCell{:}];
elseif ~isempty(obj.Prec) && isa(obj.Prec, 'aAMG')
   if obj.Prec.phys == 0, physics = "pressure"; else, physics = "displacements"; end
else
   physics = "displacements";
end

% Extract solver parameters and time step
dt = obj.generalsolver.dt;
solverType = obj.SolverType;
tol = obj.generalsolver.simparams.relTol;
maxit = obj.params.maxit;
if isfield(obj.params, 'restart') && ~isempty(obj.params.restart)
   restart = obj.params.restart;
else
   restart = 100;
end

% Persist system matrix, RHS, and full execution metadata to disk
save(filename, 'A', 'b', 'coordinates', 'ruizFlag', 'preconditioner', ...
   'physics', 'solverType', 'restart', 'tol', 'maxit', 'dt', 'R', 'timeRelaxation', 'props');
gresLog().log(3, 'Saved linear system solve %d to %s\n', currentSolve, filename);
end

% Extract spatial grid coordinates from domain objects across single or multiple domains
function coordinates = getCoordinates(generalsolver)
if numel(generalsolver.domains) == 1
   coordinates = generalsolver.domains(1).grid.coordinates;
else
   coordinates = arrayfun(@(d) d.grid.coordinates, generalsolver.domains, 'UniformOutput', false);
end
end
