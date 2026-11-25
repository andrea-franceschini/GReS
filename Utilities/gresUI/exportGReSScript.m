function scriptPath = exportGReSScript(config, targetDir, problemXmlName)
%EXPORTGRESSCRIPT Generate a standalone MATLAB script to run a configured GReS simulation.
%
%   scriptPath = exportGReSScript(config, targetDir, problemXmlName)
%
%   INPUTS:
%     config         : struct returned by gresUI.collectConfiguration()
%     targetDir      : folder where runSimulation.m is written
%     problemXmlName : filename of the master problem XML (e.g. 'problem.xml')
%
%   OUTPUT:
%     scriptPath     : path to the generated runSimulation.m script

if ~isfolder(targetDir)
  mkdir(targetDir);
end
if nargin < 3 || strlength(string(problemXmlName)) == 0
  problemXmlName = "problem.xml";
else
  [~, nm, ext] = fileparts(char(problemXmlName));
  problemXmlName = string([nm ext]);
end

lines = string.empty;
lines(end+1) = "% Auto-generated GReS simulation script";
lines(end+1) = "";
lines(end+1) = sprintf("params = readInput('%s');", problemXmlName);
lines(end+1) = "simparams = SimulationParameters(params.SimulationParameters);";

domainsCfg = config.domains;
hasGeom = ~isempty(domainsCfg) && all(arrayfun( ...
  @(d) isfield(d, 'geometry') && isstruct(d.geometry) && isfield(d.geometry, 'type'), domainsCfg));
hasInterf = isfield(config, 'interfaces') && ~isempty(config.interfaces);

schemeName = string(config.solutionScheme.name);
schemeArgs = ["'simulationparameters', simparams", "'domains', domains"];

if hasGeom
  if hasInterf
    lines(end+1) = "[domains, interfaces] = buildModel(params);";
    schemeArgs(end+1) = "'interface', interfaces";
  else
    lines(end+1) = "domains = buildModel(params);";
  end
else
  lines(end+1) = "% Define mesh grid and optional Materials/Boundaries before running:";
  lines(end+1) = "% grid = structuredMesh(NX, NY, NZ, LX, LY, LZ);";
  lines(end+1) = "domains = Discretizer('Grid', grid);";
  lines(end+1) = "domains.addPhysicsSolvers(params.Solver);";
end

hasOut = isfield(config, 'output') && isstruct(config.output) && ...
  ((isfield(config.output, 'params') && ~isempty(fieldnames(config.output.params))) || ...
   (isfield(config.output, 'restart') && ~isempty(fieldnames(config.output.restart))));

if hasOut
  lines(end+1) = "output = OutState(params.Output);";
  schemeArgs(end+1) = "'output', output";
end

if isfield(config.solutionScheme, 'params') && isstruct(config.solutionScheme.params)
  fn = fieldnames(config.solutionScheme.params);
  for i = 1:numel(fn)
    val = config.solutionScheme.params.(fn{i});
    schemeArgs(end+1) = sprintf("'%s', %s", fn{i}, formatMatlabLiteral(val)); %#ok<AGROW>
  end
end

lines(end+1) = "";
lines(end+1) = sprintf("solver = %s(%s);", schemeName, strjoin(schemeArgs, ", "));
lines(end+1) = "solver.simulationLoop();";

scriptPath = string(fullfile(char(targetDir), 'runSimulation.m'));
fid = fopen(char(scriptPath), 'w', 'n', 'UTF-8');
if fid == -1
  error('exportGReSScript:FileOpenFailed', 'Could not write "%s".', scriptPath);
end
cleanupObj = onCleanup(@() fclose(fid));
fwrite(fid, char(strjoin(lines, newline) + newline), 'char');
end


function lit = formatMatlabLiteral(val)
if islogical(val)
  lit = string(double(val));
elseif isnumeric(val)
  if isscalar(val)
    lit = string(sprintf('%.15g', val));
  else
    parts = arrayfun(@(x) sprintf('%.15g', x), val(:)', 'UniformOutput', false);
    lit = "[" + strjoin(string(parts), ", ") + "]";
  end
else
  lit = """" + string(val) + """";
end
end
