function schema = scanGReSParameters()
%SCANGREsPARAMETERS Automatically extract GReS configurable parameters from the codebase.
%
%   schema = scanGReSParameters()
%
%   Inspects SimulationParameters, PhysicsSolver subclasses, SolutionScheme
%   subclasses, InterfaceSolver/MortarQuadrature subclasses, Chronos/LinearSolver
%   options, Grid.create, OutState, Materials/ConstitutiveLaws/VanGenuchten,
%   and Boundaries directly from source files without manual schema tables.
%
%   OUTPUT:
%     schema : struct with fields:
%              .simParams, .physicsSolvers, .linearSolver, .solutionSchemes,
%              .geometry, .outputs, .materials, .boundaryConditions, .interfaces

rootDir = gres_root();

schema = struct();
schema.simParams = scanSimulationParameters(rootDir);
schema.physicsSolvers = scanPhysicsSolvers(rootDir);
schema.linearSolver = scanLinearSolver(rootDir);
schema.solutionSchemes = scanSolutionSchemes(rootDir);
schema.geometry = scanGeometry(rootDir);
schema.outputs = scanOutputs(rootDir);
schema.materials = scanMaterials(rootDir);
schema.boundaryConditions = scanBoundaryConditions(rootDir);
schema.interfaces = scanInterfaces(rootDir);

end


function params = scanSimulationParameters(rootDir)
% Parse default = struct(...) inside readSimulationParameters.
filePath = fullfile(rootDir, 'Code', 'simparam', 'SimulationParameters.m');
lines = readLines(filePath);
methodBody = extractMethodBody(lines, 'readSimulationParameters');
[params, ~] = extractReadInputParams(methodBody, 'SimulationParameters');

% LinearSolver has its own dedicated GUI section
keep = ([params.name] ~= "LinearSolver");
params = params(keep);
end


function solvers = scanPhysicsSolvers(rootDir)
% Discover concrete PhysicsSolver subclasses and extract registerSolver parameters.
solverRoot = fullfile(rootDir, 'Code', 'discretizer', 'PhysicsSolvers');
classFiles = findClassFiles(solverRoot);

solvers = struct('name', {}, 'fields', {}, 'isCoupled', {}, ...
  'hasCustomSchema', {}, 'params', {}, 'mechanicsChoices', {}, 'flowChoices', {});

allConcreteNames = string.empty;
metaMap = containers.Map();
fileMap = containers.Map();

for i = 1:numel(classFiles)
  [~, clsName] = fileparts(classFiles{i});
  mc = meta.class.fromName(clsName);
  if isempty(mc) || ~inheritsFrom(mc, 'PhysicsSolver') || mc.Abstract
    continue;
  end
  allConcreteNames(end+1) = string(clsName); %#ok<AGROW>
  metaMap(clsName) = mc;
  fileMap(clsName) = classFiles{i};
end

allConcreteNames = sort(allConcreteNames);

% Identify concrete mechanics and flow subclasses (including base concrete classes)
mechChoices = string.empty;
flowChoices = string.empty;
for i = 1:numel(allConcreteNames)
  cName = char(allConcreteNames(i));
  mc = metaMap(cName);
  if strcmp(cName, 'Poromechanics') || inheritsFrom(mc, 'Poromechanics')
    mechChoices(end+1) = string(cName); %#ok<AGROW>
  elseif strcmp(cName, 'SinglePhaseFlow') || inheritsFrom(mc, 'SinglePhaseFlow')
    flowChoices(end+1) = string(cName); %#ok<AGROW>
  end
end

for i = 1:numel(allConcreteNames)
  cName = char(allConcreteNames(i));
  mc = metaMap(cName);

  fields = string.empty;
  if hasStaticMethod(mc, 'getField')
    fields = string(feval([cName '.getField']));
  end

  [params, hasCustom] = collectSolverParams(mc, fileMap);

  isCoupled = ~isempty(params) && ...
    any([params.name] == "mechanicsSolver") && ...
    any([params.name] == "flowSolver");

  entry = struct();
  entry.name = string(cName);
  entry.fields = reshape(fields, 1, []);
  entry.isCoupled = isCoupled;
  entry.hasCustomSchema = hasCustom && isempty(params);
  entry.params = params;
  entry.mechanicsChoices = mechChoices;
  entry.flowChoices = flowChoices;

  % Coupled solvers expect child <Mechanics> and <Flow> tags, not attributes
  if isCoupled
    keep = ~ismember([entry.params.name], ["Mechanics", "Flow"]);
    entry.params = entry.params(keep);
    for pIdx = 1:numel(entry.params)
      if entry.params(pIdx).name == "mechanicsSolver"
        entry.params(pIdx).type = "enum";
        entry.params(pIdx).choices = mechChoices;
        entry.params(pIdx).allowedText = "one of [" + strjoin(mechChoices, ", ") + "]";
      elseif entry.params(pIdx).name == "flowSolver"
        entry.params(pIdx).type = "enum";
        entry.params(pIdx).choices = flowChoices;
        entry.params(pIdx).allowedText = "one of [" + strjoin(flowChoices, ", ") + "]";
      end
    end
  end

  solvers(end+1) = entry; %#ok<AGROW>
end
end


function [params, hasCustom] = collectSolverParams(mc, fileMap)
% Walk inheritance hierarchy (base -> derived) and extract registerSolver parameters.
chain = getInheritanceChain(mc);
params = emptyParamArray();
hasCustom = false;

for i = 1:numel(chain)
  clsName = chain{i};
  if strcmp(clsName, 'PhysicsSolver') || strcmp(clsName, 'handle')
    continue;
  end

  clsPath = which(clsName);
  if isempty(clsPath) && isKey(fileMap, clsName)
    clsPath = fileMap(clsName);
  end
  if isempty(clsPath) || ~isfile(clsPath)
    continue;
  end

  lines = readLines(clsPath);
  regBody = extractMethodBody(lines, 'registerSolver');
  if isempty(regBody)
    continue;
  end

  superOverrides = extractSuperclassOverrides(regBody, 'registerSolver');
  params = applyOverrides(params, superOverrides);

  helperNames = extractVararginHelperCalls(regBody, clsName);
  combinedBody = regBody;
  for h = 1:numel(helperNames)
    hBody = extractMethodBody(lines, helperNames{h});
    if ~isempty(hBody)
      combinedBody = [combinedBody; hBody]; %#ok<AGROW>
    end
  end

  delegatedSolvers = extractDelegatedSolvers(regBody);
  for d = 1:numel(delegatedSolvers)
    dMc = meta.class.fromName(delegatedSolvers{d});
    if ~isempty(dMc)
      [dParams, ~] = collectSolverParams(dMc, fileMap);
      params = mergeParamArrays(params, dParams);
    end
  end

  [localParams, foundReadInput] = extractReadInputParams(combinedBody, clsName);
  params = mergeParamArrays(params, localParams);

  if ~foundReadInput && any(contains(regBody, 'checkInput'))
    hasCustom = true;
  end
end
end


function schemes = scanSolutionSchemes(rootDir)
% Discover concrete SolutionScheme subclasses and extract setSolutionScheme options.
schemeRoot = fullfile(rootDir, 'Code', 'problemsolver', 'SolutionSchemes');
classFiles = findClassFiles(schemeRoot);
schemes = struct('name', {}, 'params', {});

skipObjectKeys = ["simulationparameters", "output", "domains", "interface"];
allNames = string.empty;
fileMap = containers.Map();

for i = 1:numel(classFiles)
  [~, clsName] = fileparts(classFiles{i});
  mc = meta.class.fromName(clsName);
  if isempty(mc) || ~inheritsFrom(mc, 'SolutionScheme') || mc.Abstract
    continue;
  end
  allNames(end+1) = string(clsName); %#ok<AGROW>
  fileMap(clsName) = classFiles{i};
end
allNames = sort(allNames);

for i = 1:numel(allNames)
  cName = char(allNames(i));
  lines = readLines(fileMap(cName));
  body = extractMethodBody(lines, 'setSolutionScheme');
  [params, ~] = extractReadInputParams(body, cName);

  if ~isempty(params)
    keep = ~ismember([params.name], skipObjectKeys);
    params = params(keep);
  end

  schemes(end+1) = struct('name', string(cName), 'params', params); %#ok<AGROW>
end
end


function geomSchema = scanGeometry(rootDir)
% Parse Grid.create modes and parameters from Code/grid/@Grid/Grid.m.
gridPath = fullfile(rootDir, 'Code', 'grid', '@Grid', 'Grid.m');
lines = readLines(gridPath);
body = extractMethodBody(lines, 'create');
[cleanLines, ~] = collapseContinuations(body);

% Top-level modes inside Grid.create
geomSchema = struct('name', {}, 'params', {});
modes = ["StructuredMesh", "BlockStructuredMesh", "meshFile"];

for m = 1:numel(modes)
  mName = modes(m);
  if mName == "meshFile"
    p = emptyParam();
    p.name = "meshFile";
    p.type = "file";
    p.required = true;
    p.defaultValue = "";
    p.defaultText = "none (compulsory)";
    p.allowedText = "path to mesh file (.msh, .vtk)";
    p.description = "External mesh file path imported by Grid.create";
    p.source = "Grid.create";
    geomSchema(end+1) = struct('name', mName, 'params', p); %#ok<AGROW>
    continue;
  end

  % Extract the branch lines for StructuredMesh / BlockStructuredMesh
  branchLines = extractBranchLines(cleanLines, char(mName));
  [params, ~] = extractReadInputParams(branchLines, char(mName));
  for k = 1:numel(params)
    if ismember(params(k).name, ["LX", "LY", "LZ"])
      params(k).type = "vector";
      params(k).allowedText = "2-element extent [min,max] (e.g. 0,1)";
    end
  end
  geomSchema(end+1) = struct('name', mName, 'params', params); %#ok<AGROW>
end
end


function outSchema = scanOutputs(rootDir)
% Parse OutState.setOutState parameters and nested restart defaults.
outPath = fullfile(rootDir, 'Code', 'output', 'OutState.m');
lines = readLines(outPath);
body = extractMethodBody(lines, 'setOutState');
[cleanLines, ~] = collapseContinuations(body);

[allParams, ~] = extractStructVarParams(cleanLines, 'default', 'OutState');
keep = ([allParams.name] ~= "restart");
params = allParams(keep);

% Extract validFormats for vtkFormat
bodyText = strjoin(cleanLines, newline);
fmtTok = regexp(bodyText, 'validFormats\s*=\s*\[([^\]]+)\]', 'tokens', 'once');
vtkChoices = ["auto", "ascii", "binary"];
if ~isempty(fmtTok)
  toks = regexp(fmtTok{1}, '[''"]([^''"]+)[''"]', 'tokens');
  if ~isempty(toks)
    vtkChoices = string(cellfun(@(c) c{1}, toks, 'UniformOutput', false));
  end
end

for k = 1:numel(params)
  switch params(k).name
    case "vtkFormat"
      params(k).type = "enum";
      params(k).choices = vtkChoices;
      params(k).allowedText = "one of [" + strjoin(vtkChoices, ", ") + "]";
    case {"outputFile", "matFileName"}
      params(k).type = "string";
      params(k).allowedText = "file/folder path prefix";
    case "printTimes"
      params(k).type = "vector";
      params(k).allowedText = "comma-separated times (e.g. 0,5,10)";
    case {"saveHistory", "solvePrintTimes"}
      params(k).type = "flag01";
      params(k).allowedText = "0 or 1 (boolean flag)";
      if isempty(params(k).defaultValue)
        params(k).defaultValue = 1;
        params(k).defaultText = "1";
      end
  end
end

[restartParams, ~] = extractStructVarParams(cleanLines, 'restartDefault', 'OutState.restart');
for k = 1:numel(restartParams)
  if restartParams(k).name == "filename"
    restartParams(k).type = "string";
    restartParams(k).allowedText = "restart file path";
  end
end

outSchema = struct('params', params, 'restartParams', restartParams);
end


function matSchema = scanMaterials(rootDir)
% Extract Fluid, Solid, PorousRock, VanGenuchten, and ConstitutiveLaw schemas.
matRoot = fullfile(rootDir, 'Code', 'materials');

% Fluid
fluidLines = readLines(fullfile(matRoot, 'Fluid.m'));
fluidBody = extractMethodBody(fluidLines, 'readMaterialParameters');
[fluidParams, ~] = extractReadInputParams(fluidBody, 'Fluid');

% Solid header
matLines = readLines(fullfile(matRoot, 'Materials.m'));
solidBody = extractMethodBody(matLines, 'addSolid');
[solidParams, ~] = extractReadInputParams(solidBody, 'Solid');
solidParams = solidParams(~ismember([solidParams.name], ["Constitutive", "PorousRock"]));
for k = 1:numel(solidParams)
  if solidParams(k).name == "cellTags"
    solidParams(k).type = "vector";
    solidParams(k).allowedText = "comma-separated cell tags (e.g. 1,2)";
  elseif solidParams(k).name == "name"
    solidParams(k).type = "string";
    solidParams(k).defaultValue = "mat_1";
    solidParams(k).defaultText = "mat_1";
    solidParams(k).allowedText = "material name";
  end
end

% PorousRock
prLines = readLines(fullfile(matRoot, 'PorousRock.m'));
prBody = extractMethodBody(prLines, 'readMaterialParameters');
[porousRockParams, ~] = extractReadInputParams(prBody, 'PorousRock');
porousRockParams = porousRockParams([porousRockParams.name] ~= "Curves");
for k = 1:numel(porousRockParams)
  if porousRockParams(k).name == "permeability"
    porousRockParams(k).type = "vector";
    porousRockParams(k).allowedText = "1, 3, or 6 numeric values (e.g. 1e-12)";
  end
end

% Constitutive laws under Code/materials/ConstitutiveLaws/
lawFiles = findClassFiles(fullfile(matRoot, 'ConstitutiveLaws'));
laws = struct('name', {}, 'params', {});
lawNames = string.empty;
lawFileMap = containers.Map();
lawSuperMap = containers.Map();
for i = 1:numel(lawFiles)
  [~, cName] = fileparts(lawFiles{i});
  lines = readLines(lawFiles{i});
  [isAbs, supers] = parseClassHeader(lines);
  if isAbs
    continue;
  end
  lawSuperMap(cName) = supers;
  txt = strjoin(lines, newline);
  hasDefaultReadInput = ~isempty(regexp(txt, 'readInput\s*\(\s*default\>', 'once'));
  if hasDefaultReadInput && (any(ismember(supers, {'ConstitutiveLaw', 'Elastic'})) || strcmp(cName, 'SedimentMaterial'))
    lawNames(end+1) = string(cName); %#ok<AGROW>
    lawFileMap(cName) = lawFiles{i};
  end
end
lawNames = sort(lawNames);

for i = 1:numel(lawNames)
  cName = char(lawNames(i));
  lawParams = collectConstitutiveLawParams(cName, lawFileMap, lawSuperMap);
  if strcmp(cName, 'SedimentMaterial')
    for k = 1:numel(lawParams)
      if lawParams(k).name == "conductivity"
        lawParams(k).type = "vector";
        lawParams(k).allowedText = "1 or 3 numeric values";
      end
    end
  end
  laws(end+1) = struct('name', string(cName), 'params', lawParams); %#ok<AGROW>
end

% VanGenuchten capillary curves
vgLines = readLines(fullfile(matRoot, 'VanGenuchten.m'));
vgBody = extractMethodBody(vgLines, 'readInputParameters');
[vgClean, ~] = collapseContinuations(vgBody);
vgText = strjoin(vgClean, newline);
caseMatches = regexp(vgText, 'case\s+[''"]([^''"]+)[''"]', 'tokens');
vgCurves = struct('type', {}, 'params', {});
for i = 1:numel(caseMatches)
  cType = string(caseMatches{i}{1});
  branchLines = extractCaseBranchLines(vgClean, char(cType));
  [cParams, ~] = extractReadInputParams(branchLines, ['VanGenuchten.' char(cType)]);
  for k = 1:numel(cParams)
    if ismember(cParams(k).name, ["soilName", "capillaryCurvePath", "relativePermeabilityPath"])
      cParams(k).type = "string";
      cParams(k).allowedText = "string / file path";
    end
  end
  vgCurves(end+1) = struct('type', cType, 'params', cParams); %#ok<AGROW>
end

matSchema = struct( ...
  'fluidParams', fluidParams, ...
  'solidParams', solidParams, ...
  'constitutiveLaws', laws, ...
  'porousRockParams', porousRockParams, ...
  'vanGenuchten', vgCurves);
end


function params = collectConstitutiveLawParams(clsName, fileMap, superMap)
% Walk file-based inheritance chain and collect readInput defaults from constructor and readMaterialParameters.
chain = getFileInheritanceChain(clsName, superMap);
params = emptyParamArray();
for i = 1:numel(chain)
  cName = chain{i};
  if ismember(cName, {'ConstitutiveLaw', 'matlab.mixin.Copyable', 'handle'})
    continue;
  end
  if isKey(fileMap, cName)
    clsPath = fileMap(cName);
  else
    clsPath = which(cName);
  end
  if isempty(clsPath) || ~isfile(clsPath)
    continue;
  end
  lines = readLines(clsPath);
  ctorBody = extractMethodBody(lines, cName);
  readBody = extractMethodBody(lines, 'readMaterialParameters');
  combined = [ctorBody; readBody];
  [localParams, ~] = extractReadInputParams(combined, cName);
  params = mergeParamArrays(params, localParams);
end
end


function [isAbs, supers] = parseClassHeader(lines)
isAbs = false;
supers = {};
for i = 1:numel(lines)
  code = strtrim(stripComment(lines{i}));
  if isempty(code)
    continue;
  end
  if startsWith(code, 'classdef')
    isAbs = contains(code, 'Abstract');
    tok = regexp(code, '<\s*(.+)$', 'tokens', 'once');
    if ~isempty(tok)
      parts = regexp(tok{1}, '&', 'split');
      supers = cellfun(@strtrim, parts, 'UniformOutput', false);
    end
    return;
  end
end
end


function chain = getFileInheritanceChain(clsName, superMap)
chain = {};
if isKey(superMap, clsName)
  supers = superMap(clsName);
  for i = 1:numel(supers)
    sub = getFileInheritanceChain(supers{i}, superMap);
    chain = [chain, sub]; %#ok<AGROW>
  end
end
chain{end+1} = clsName;
chain = unique(chain, 'stable');
end


function bcSchema = scanBoundaryConditions(rootDir)
% Extract addBC and addBCEvent parameters plus enum choices from entityField/BCtype/readEntitySet.
addBcPath = fullfile(rootDir, 'Code', 'boundaries', '@Boundaries', 'addBC.m');
bcLines = readLines(addBcPath);
[bcParams, ~] = extractReadInputParams(bcLines, 'Boundaries.addBC');
bcParams = bcParams([bcParams.name] ~= "options");

% Enum choices directly from MATLAB enumerations and readEntitySet
fieldEnums = string(enumeration('entityField'))';
fieldEnums = fieldEnums(fieldEnums ~= "face");
typeEnums = string(enumeration('BCtype'))';

entSetLines = readLines(fullfile(rootDir, 'Code', 'boundaries', '@BoundaryEntities', 'readEntitySet.m'));
entText = strjoin(entSetLines, newline);
entTypes = ["tags", "bclist", "box", "bclistfile"];
if contains(entText, 'case "bclistfile"')
  % Keep canonical choices in readable order
  entTypes = ["tags", "bclist", "box", "bclistfile"];
end

for k = 1:numel(bcParams)
  switch bcParams(k).name
    case {"name", "variable"}
      bcParams(k).type = "string";
      bcParams(k).allowedText = "string identifier";
    case "field"
      bcParams(k).type = "enum";
      bcParams(k).choices = fieldEnums;
      bcParams(k).defaultValue = fieldEnums(1);
      bcParams(k).defaultText = "none (compulsory)";
      bcParams(k).allowedText = "one of [" + strjoin(fieldEnums, ", ") + "]";
    case "entityListType"
      bcParams(k).type = "enum";
      bcParams(k).choices = entTypes;
      bcParams(k).defaultValue = entTypes(1);
      bcParams(k).defaultText = "none (compulsory)";
      bcParams(k).allowedText = "one of [" + strjoin(entTypes, ", ") + "]";
    case "entityList"
      bcParams(k).type = "string";
      bcParams(k).allowedText = "comma-separated IDs/tags/box (e.g. 1,2) or file path";
      bcParams(k).helpText = strjoin([ ...
        "Specifies which mesh entities receive this boundary condition (interpreted using 'entityListType' and 'field').", ...
        "", ...
        "- tags: surface/cell tag IDs (e.g. 1 or 1,2). With GReS StructuredMesh (field = surface):", ...
        "    1 = Bottom (z min),   2 = Top (z max)", ...
        "    3 = South (y min),    4 = North (y max)", ...
        "    5 = West (x min),     6 = East (x max)", ...
        "    (StructuredMesh cells use tag 1; in octree-refined BlockStructuredMesh with non-conforming faces, tags 1-2 are internal faces and external faces shift to 3-8.)", ...
        "- bclist: explicit 1-based entity indices (node/cell/surface IDs).", ...
        "- box: 6 coordinates [xmin, xmax, ymin, ymax, zmin, zmax] selecting entities strictly inside the box.", ...
        "- bclistfile: path to a text file listing entity indices."], newline);
    case "type"
      bcParams(k).type = "editableEnum";
      bcParams(k).choices = typeEnums;
      bcParams(k).defaultValue = typeEnums(1);
      bcParams(k).defaultText = "none (compulsory)";
      bcParams(k).allowedText = "one of [" + strjoin(typeEnums, ", ") + "] or custom type";
    case "components"
      bcParams(k).type = "string";
      bcParams(k).allowedText = "e.g. x, y, z or x,y";
    case "essential"
      bcParams(k).type = "flag01";
      bcParams(k).allowedText = "0 or 1 (boolean flag)";
  end
end

evPath = fullfile(rootDir, 'Code', 'boundaries', '@BoundaryEntities', 'addBCEvent.m');
evLines = readLines(evPath);
[eventParams, ~] = extractReadInputParams(evLines, 'BoundaryEntities.addBCEvent');

bcSchema = struct('bcParams', bcParams, 'eventParams', eventParams);
end


function interfSchema = scanInterfaces(rootDir)
% Discover InterfaceSolver and MortarQuadrature subclasses under Code/multidomain.
multRoot = fullfile(rootDir, 'Code', 'multidomain');

% Common Mortar interface parameters from InterfaceSolver & BBTree
basePath = fullfile(multRoot, 'InterfaceSolver.m');
baseLines = readLines(basePath);
ctorBody = extractMethodBody(baseLines, 'InterfaceSolver');
mortarBody = extractMethodBody(baseLines, 'setMortarInterface');
procBody = extractMethodBody(baseLines, 'processMortarGrid');
varBody = extractMethodBody(baseLines, 'setCoupledVariables');

[ctorParams, ~] = extractReadInputParams(ctorBody, 'InterfaceSolver');
[mortarParams, ~] = extractReadInputParams(mortarBody, 'InterfaceSolver.setMortarInterface');
if ~isempty(mortarParams)
  mortarParams = mortarParams([mortarParams.name] ~= "Quadrature");
end
[varParams, ~] = extractReadInputParams(varBody, 'InterfaceSolver.setCoupledVariables');

bbLines = readLines(fullfile(multRoot, 'Mortar', 'BBTree.m'));
bbCtor = extractMethodBody(bbLines, 'BBTree');
[bbParams, ~] = extractReadInputParams(bbCtor, 'BBTree');
if ~isempty(bbParams)
  bbParams = bbParams([bbParams.name] == "scale");
end

rawCommon = mergeParamArrays(ctorParams, mortarParams);
rawCommon = mergeParamArrays(rawCommon, varParams);
rawCommon = mergeParamArrays(rawCommon, bbParams);

% Extract multiplierType choices from processMortarGrid switch block
multChoices = ["P0", "standard", "dual"];
procText = strjoin(procBody, newline);
caseToks = regexp(procText, 'case\s+\{?([^\n\r\}]+)\}?', 'tokens');
foundMult = string.empty;
for k = 1:numel(caseToks)
  sTok = regexp(caseToks{k}{1}, '[''"]([^''"]+)[''"]', 'tokens');
  for m = 1:numel(sTok)
    foundMult(end+1) = string(sTok{m}{1}); %#ok<AGROW>
  end
end
if ~isempty(foundMult)
  multChoices = unique(foundMult, 'stable');
end

commonOrder = ["masterDomain", "slaveDomain", "masterSurface", "slaveSurface", ...
  "multiplierType", "variable", "scale"];
commonParams = emptyParamArray();
for k = 1:numel(commonOrder)
  idx = find([rawCommon.name] == commonOrder(k), 1);
  if isempty(idx)
    continue;
  end
  p = rawCommon(idx);
  switch p.name
    case "masterDomain"
      p.type = "domainRef";
      p.required = true;
      p.defaultValue = "2";
      p.defaultText = "2";
      p.allowedText = "domain index (1..nDom)";
      p.description = "Master domain index for mortar coupling";
    case "slaveDomain"
      p.type = "domainRef";
      p.required = true;
      p.defaultValue = "1";
      p.defaultText = "1";
      p.allowedText = "domain index (1..nDom)";
      p.description = "Slave domain index for mortar coupling";
    case {"masterSurface", "slaveSurface"}
      p.type = "vector";
      p.required = false;
      p.defaultValue = double.empty;
      p.defaultText = "all boundary surfaces";
      p.allowedText = "surface tag ID(s) (e.g. 5 or 5,6)";
    case "multiplierType"
      p.type = "enum";
      p.choices = multChoices;
      p.defaultValue = "P0";
      p.defaultText = "P0";
      p.allowedText = "one of [" + strjoin(multChoices, ", ") + "]";
    case "variable"
      p.type = "editableEnum";
      p.required = false;
      p.choices = ["displacements", "pressure"];
      p.defaultValue = "";
      p.defaultText = "auto (shared domain variable)";
      p.allowedText = "coupled variable name or leave unchecked for auto";
    case "scale"
      p.type = "double";
      p.required = false;
      p.allowedText = "positive scalar (BBTree search inflation)";
  end
  commonParams(end+1) = p; %#ok<AGROW>
end

% Quadrature subclasses under Code/multidomain/Mortar/Quadrature
quadRoot = fullfile(multRoot, 'Mortar', 'Quadrature');
mqLines = readLines(fullfile(quadRoot, 'MortarQuadrature.m'));
mqCtor = extractMethodBody(mqLines, 'MortarQuadrature');
[baseQuadParams, ~] = extractReadInputParams(mqCtor, 'MortarQuadrature');

quadFiles = findClassFiles(quadRoot);
quadratures = struct('name', {}, 'params', {});
quadNames = string.empty;
quadFileMap = containers.Map();
for i = 1:numel(quadFiles)
  [~, qName] = fileparts(quadFiles{i});
  mc = meta.class.fromName(qName);
  if isempty(mc) || ~inheritsFrom(mc, 'MortarQuadrature') || mc.Abstract
    continue;
  end
  quadNames(end+1) = string(qName); %#ok<AGROW>
  quadFileMap(qName) = quadFiles{i};
end
quadNames = sort(quadNames);
% Put SegmentBasedQuadrature first as GReS default
segIdx = find(quadNames == "SegmentBasedQuadrature", 1);
if ~isempty(segIdx) && segIdx > 1
  quadNames = ["SegmentBasedQuadrature", quadNames(quadNames ~= "SegmentBasedQuadrature")];
end

for i = 1:numel(quadNames)
  qName = char(quadNames(i));
  qLines = readLines(quadFileMap(qName));
  qCtor = extractMethodBody(qLines, qName);
  [subQParams, ~] = extractReadInputParams(qCtor, qName);
  qParams = mergeParamArrays(baseQuadParams, subQParams);
  if strcmp(qName, 'RBFquadrature')
    rbfIdx = find([qParams.name] == "RBFtype", 1);
    if ~isempty(rbfIdx)
      rbfText = strjoin(qLines, newline);
      cToks = regexp(rbfText, 'case\s+[''"]([^''"]+)[''"]', 'tokens');
      rbfChoices = ["gauss", "wendland"];
      if ~isempty(cToks)
        rbfChoices = unique(string(cellfun(@(c) c{1}, cToks, 'UniformOutput', false)), 'stable');
      end
      qParams(rbfIdx).type = "enum";
      qParams(rbfIdx).choices = rbfChoices;
      qParams(rbfIdx).allowedText = "one of [" + strjoin(rbfChoices, ", ") + "]";
    end
  end
  quadratures(end+1) = struct('name', string(qName), 'params', qParams); %#ok<AGROW>
end

% ActiveSet defaults (shared helper for contact interfaces)
asPath = fullfile(multRoot, 'interfaceModels', 'Contact', 'initializeActiveSet.m');
asLines = readLines(asPath);
[asClean, ~] = collapseContinuations(asLines);
asText = strjoin(asClean, newline);
[asParams, ~] = extractStructVarParams(asClean, 'params', 'ActiveSet');
if isempty(asParams)
  tokAs = regexp(asText, 'readInput\s*\(\s*struct\s*\(([^)]+)\)\s*,\s*activeSetParams', 'tokens', 'once');
  if ~isempty(tokAs)
    asParams = parseStructArgs(tokAs{1}, containers.Map(), 'ActiveSet', '');
  end
end
if ~isempty(asParams)
  asParams = asParams([asParams.name] ~= "Tolerances");
end
for k = 1:numel(asParams)
  if asParams(k).name == "resetActiveSet"
    asParams(k).type = "flag01";
    asParams(k).allowedText = "0 or 1 (boolean flag)";
  elseif asParams(k).name == "forceStickBoundary"
    asParams(k).type = "string";
    asParams(k).defaultValue = "";
    asParams(k).defaultText = "none (optional)";
    asParams(k).allowedText = "surface tag, axis, or empty";
  end
end
[asTolParams, ~] = extractStructVarParams(asClean, 'default', 'ActiveSet.Tolerances');
if ~isempty(asTolParams)
  excludeTol = "Tolerances";
  if ~isempty(asParams)
    excludeTol = [excludeTol, [asParams.name]];
  end
  asTolParams = asTolParams(~ismember([asTolParams.name], excludeTol));
end

% Concrete InterfaceSolver subclasses
interfFiles = findClassFiles(multRoot);
interfNames = string.empty;
interfMeta = containers.Map();
interfFileMap = containers.Map();
for i = 1:numel(interfFiles)
  [~, iName] = fileparts(interfFiles{i});
  mc = meta.class.fromName(iName);
  if isempty(mc) || ~inheritsFrom(mc, 'InterfaceSolver') || mc.Abstract
    continue;
  end
  interfNames(end+1) = string(iName); %#ok<AGROW>
  interfMeta(iName) = mc;
  interfFileMap(iName) = interfFiles{i};
end
interfNames = sort(interfNames);

models = struct( ...
  'name', {}, ...
  'params', {}, ...
  'supportsQuadrature', {}, ...
  'supportsMultiplierType', {}, ...
  'multiplierChoices', {}, ...
  'requiredField', {}, ...
  'coulombParams', {}, ...
  'activeSetParams', {}, ...
  'activeSetTolParams', {});

for i = 1:numel(interfNames)
  iName = char(interfNames(i));
  mc = interfMeta(iName);
  chain = getInheritanceChain(mc);

  mParams = emptyParamArray();
  hasCoulomb = false;
  hasActiveSet = false;
  coulombParams = emptyParamArray();

  for c = 1:numel(chain)
    cName = chain{c};
    if ismember(cName, {'InterfaceSolver', 'handle'})
      continue;
    end
    cPath = which(cName);
    if isempty(cPath) && isKey(interfFileMap, cName)
      cPath = interfFileMap(cName);
    end
    if isempty(cPath) || ~isfile(cPath)
      continue;
    end
    cLines = readLines(cPath);
    regBody = extractMethodBody(cLines, 'registerInterface');
    if isempty(regBody)
      continue;
    end
    [localParams, ~] = extractReadInputParams(regBody, cName);
    if ~isempty(localParams)
      localParams = localParams(~ismember([localParams.name], ["Coulomb", "ActiveSet", "Quadrature"]));
      mParams = mergeParamArrays(mParams, localParams);
    end

    regText = strjoin(regBody, newline);
    if contains(regText, 'input.Coulomb')
      hasCoulomb = true;
      [cleanReg, ~] = collapseContinuations(regBody);
      cleanRegText = strjoin(cleanReg, newline);
      cTok = regexp(cleanRegText, 'readInput\s*\(\s*struct\s*\(([^)]+)\)\s*,\s*input\.Coulomb', 'tokens', 'once');
      if ~isempty(cTok)
        coulombParams = parseStructArgs(cTok{1}, containers.Map(), [cName '.Coulomb'], '');
      end
    end
    if contains(regText, 'initializeActiveSet')
      hasActiveSet = true;
    end
  end

  if hasCoulomb && ~isempty(coulombParams) && ~isempty(mParams)
    mParams = mParams(~ismember([mParams.name], [coulombParams.name]));
  end

  for k = 1:numel(mParams)
    if ismember(mParams(k).name, ["forceStick", "isSmooth"])
      mParams(k).type = "flag01";
      mParams(k).allowedText = "0 or 1 (boolean flag)";
    end
  end

  srcText = fileread(interfFileMap(iName));
  supportsQuad = isempty(regexp(srcText, 'params\.Quadrature\s*=', 'once'));
  supportsMult = isempty(regexp(srcText, 'params\.multiplierType\s*=', 'once'));

  allowedMult = multChoices;
  if ~isempty(regexp(srcText, 'multiplierLocation\s*~=\s*entityField\.surface', 'once'))
    allowedMult = "P0";
  end

  reqField = "";
  if ~isempty(regexp(srcText, 'strcmp\([^,]+,\s*"pressure"\)', 'once'))
    reqField = "pressure";
  elseif ~isempty(regexp(srcText, 'getState\([^)]*"displacements"', 'once'))
    reqField = "displacements";
  end

  entry = struct();
  entry.name = string(iName);
  entry.params = mParams;
  entry.supportsQuadrature = supportsQuad;
  entry.supportsMultiplierType = supportsMult;
  entry.multiplierChoices = allowedMult;
  entry.requiredField = reqField;
  if hasCoulomb
    entry.coulombParams = coulombParams;
  else
    entry.coulombParams = emptyParamArray();
  end
  if hasActiveSet
    entry.activeSetParams = asParams;
    entry.activeSetTolParams = asTolParams;
  else
    entry.activeSetParams = emptyParamArray();
    entry.activeSetTolParams = emptyParamArray();
  end
  models(end+1) = entry; %#ok<AGROW>
end

interfSchema = struct( ...
  'commonParams', commonParams, ...
  'models', models, ...
  'quadratures', quadratures);
end


function linSchema = scanLinearSolver(rootDir)
% Extract top-level linSolverParams options and permitted Chronos AMG parameters.
linSchema = struct();
linSchema.topLevel = scanTopLevelLinearParams(rootDir);
linSchema.chronos = scanChronosParameters(rootDir);
end


function topParams = scanTopLevelLinearParams(rootDir)
% Discover isfield(..., '<param>') checks on linSolverParams across Code/linsolver.
linRoot = fullfile(rootDir, 'Code', 'linsolver');
mFiles = findMFiles(linRoot);

ignoreFields = {'Flow', 'Mechanics', 'amg', 'smoother', 'prolong', 'coarsen', 'tspace'};
paramMap = containers.Map();

for i = 1:numel(mFiles)
  lines = readLines(mFiles{i});
  joined = strjoin(lines, newline);
  if ~contains(joined, 'linSolverParams')
    continue;
  end

  aliases = {'linSolverParams'};
  aliasMatches = regexp(joined, '\<(\w+)\s*=\s*[\w\.]*simparams\.linSolverParams\>', 'tokens');
  for k = 1:numel(aliasMatches)
    aliases{end+1} = aliasMatches{k}{1}; %#ok<AGROW>
  end
  aliases = unique(aliases);

  [~, shortName] = fileparts(mFiles{i});

  for ln = 1:numel(lines)
    rawLine = lines{ln};
    codePart = stripComment(rawLine);

    for a = 1:numel(aliases)
      pat = ['isfield\([^,]*\<' aliases{a} '\>\s*,\s*[''"](\w+)[''"]\s*\)'];
      tok = regexp(codePart, pat, 'tokens');
      for t = 1:numel(tok)
        fName = tok{t}{1};
        if any(strcmp(fName, ignoreFields))
          continue;
        end
        if ~isKey(paramMap, fName)
          p = inferTopLevelParam(fName, lines, ln, shortName);
          paramMap(fName) = p;
        end
      end
    end
  end
end

keysList = sort(keys(paramMap));
topParams = emptyParamArray();
for i = 1:numel(keysList)
  topParams(end+1) = paramMap(keysList{i}); %#ok<AGROW>
end
end


function p = inferTopLevelParam(fName, lines, lineIdx, sourceName)
% Infer type, default value, choices, and description from surrounding lines.
p = emptyParam();
p.name = string(fName);
p.source = string(sourceName);
p.required = false;

ctxStart = max(1, lineIdx - 4);
ctxEnd = min(numel(lines), lineIdx + 12);
ctxLines = lines(ctxStart:ctxEnd);
ctxText = strjoin(ctxLines, newline);

propDefault = findPropertyDefault(lines, fName);
if isempty(propDefault) && strcmp(fName, 'tol')
  propDefault = findPropertyDefault(lines, 'nonLinTol');
end

desc = extractLeadingComment(lines, lineIdx);
if isempty(desc)
  desc = sprintf('Linear solver option (%s)', sourceName);
end
p.description = string(desc);

if contains(lower(fName), 'solves')
  p.type = "vector";
  p.defaultValue = double.empty;
  p.defaultText = "empty (optional)";
  p.allowedText = "numeric indices (e.g. 1,2,3)";
  return;
end

choices = string.empty;
if ~isempty(regexp(ctxText, ['switch\s+[\w\.]*\<' fName '\>'], 'once'))
  caseTokens = regexp(ctxText, 'case\s+([^\n]+)', 'tokens');
  for c = 1:numel(caseTokens)
    strToks = regexp(caseTokens{c}{1}, '[''"]([^''"]+)[''"]', 'tokens');
    for s = 1:numel(strToks)
      choices(end+1) = string(strToks{s}{1}); %#ok<AGROW>
    end
  end
end

if ~isempty(choices)
  p.type = "enum";
  p.choices = unique(choices, 'stable');
  if ~isempty(propDefault) && isstring(propDefault)
    p.defaultValue = propDefault;
  else
    p.defaultValue = p.choices(1);
  end
  p.defaultText = p.defaultValue;
  p.allowedText = "one of [" + strjoin(p.choices, ", ") + "]";
  return;
end

isFlag = startsWith(fName, 'use') || endsWith(fName, 'Flag') || ...
  contains(ctxText, '== 1') || contains(ctxText, 'logical(');

if isFlag
  p.type = "flag01";
  if ~isempty(propDefault) && islogical(propDefault)
    p.defaultValue = double(propDefault);
  elseif ~isempty(propDefault) && isnumeric(propDefault)
    p.defaultValue = double(propDefault ~= 0);
  else
    p.defaultValue = 0;
  end
  p.defaultText = string(p.defaultValue);
  p.allowedText = "0 or 1 (boolean flag)";
  return;
end

if ~isempty(propDefault) && (isstring(propDefault) || ischar(propDefault))
  p.type = "string";
  p.defaultValue = string(propDefault);
  p.defaultText = p.defaultValue;
  p.allowedText = "string";
  return;
end

p.type = "double";
if ~isempty(propDefault) && isnumeric(propDefault)
  p.defaultValue = propDefault;
  p.defaultText = string(propDefault);
else
  p.defaultValue = double.empty;
  p.defaultText = "none (optional)";
end
p.allowedText = "real scalar";
end


function chronos = scanChronosParameters(rootDir)
% Read Chronos XML defaults and filter out blocked fields from @aAMG/getUserInput.m.
getUserInputPath = fullfile(rootDir, 'Code', 'linsolver', 'preconditioning', '@aAMG', 'getUserInput.m');
uiLines = readLines(getUserInputPath);
uiText = strjoin(uiLines, newline);

secMatches = regexp(uiText, 'isfield\(\s*input\s*,\s*[''"](\w+)[''"]\s*\)', 'tokens');
allowedSections = string.empty;
for i = 1:numel(secMatches)
  allowedSections(end+1) = string(secMatches{i}{1}); %#ok<AGROW>
end
allowedSections = unique(allowedSections, 'stable');

blockedMap = containers.Map();
for i = 1:numel(allowedSections)
  sec = char(allowedSections(i));
  capSec = [upper(sec(1)) sec(2:end)];
  pat = ['blocked' capSec '\s*=\s*\{([^}]*)\}'];
  joinedClean = regexprep(uiText, '\.\.\.\s*\n', ' ');
  tok = regexp(joinedClean, pat, 'tokens', 'once');
  if ~isempty(tok)
    items = regexp(tok{1}, '[''"]([^''"]+)[''"]', 'tokens');
    blockedMap(sec) = cellfun(@(c) c{1}, items, 'UniformOutput', false);
  else
    blockedMap(sec) = {};
  end
end

flowXml = fullfile(rootDir, 'Code', 'linsolver', 'XML_setup', 'chronos_xml_setup_CFD.xml');
mechXml = fullfile(rootDir, 'Code', 'linsolver', 'XML_setup', 'chronos_xml_setup.xml');

chronos = struct();
chronos.Flow = buildChronosPhysicsSchema(flowXml, allowedSections, blockedMap, 'Flow');
chronos.Mechanics = buildChronosPhysicsSchema(mechXml, allowedSections, blockedMap, 'Mechanics');
end


function physSchema = buildChronosPhysicsSchema(xmlPath, allowedSections, blockedMap, physLabel)
rawData = readInput(xmlPath);
physSchema = struct('section', {}, 'params', {});

for i = 1:numel(allowedSections)
  sec = char(allowedSections(i));
  if ~isfield(rawData, sec)
    continue;
  end

  secStruct = rawData.(sec);
  fNames = fieldnames(secStruct);
  blocked = blockedMap(sec);

  params = emptyParamArray();
  for k = 1:numel(fNames)
    fn = fNames{k};
    if any(strcmpi(fn, blocked))
      continue;
    end

    val = secStruct.(fn);
    p = emptyParam();
    p.name = string(fn);
    p.required = false;
    p.source = string([physLabel '.' sec]);
    p.description = sprintf('Chronos %s %s parameter', physLabel, sec);

    if isstring(val) || ischar(val)
      p.type = "string";
      p.defaultValue = string(val);
      p.defaultText = string(val);
      p.allowedText = "string";
    else
      p.type = "double";
      p.defaultValue = double(val);
      p.defaultText = string(val);
      p.allowedText = "numeric scalar";
    end
    params(end+1) = p; %#ok<AGROW>
  end

  if ~isempty(params)
    entry = struct();
    entry.section = string(sec);
    entry.params = params;
    physSchema(end+1) = entry; %#ok<AGROW>
  end
end
end


function [params, foundReadInput] = extractReadInputParams(lines, sourceName)
% Extract key-value defaults from `default = struct(...)` and inline `readInput(struct(...), ...)`.
[cleanLines, ~] = collapseContinuations(lines);
[params, foundReadInput] = extractStructVarParams(cleanLines, 'default', sourceName);
end


function [params, foundReadInput] = extractStructVarParams(cleanLines, varName, sourceName)
% Extract key-value pairs from `<varName> = struct(...)` and inline `readInput(struct(...), ...)`.
params = emptyParamArray();
foundReadInput = false;

localVars = containers.Map();
for i = 1:numel(cleanLines)
  assignTok = regexp(cleanLines{i}, '^\s*(\w+)\s*=\s*([^;]+);?\s*$', 'tokens', 'once');
  if ~isempty(assignTok) && ~strcmp(assignTok{1}, varName)
    localVars(assignTok{1}) = strtrim(assignTok{2});
  end
end

for i = 1:numel(cleanLines)
  line = cleanLines{i};
  if contains(line, 'readInput(')
    foundReadInput = true;
  end

  tokDef = regexp(line, ['\<' varName '\s*=\s*struct\s*\((.*)\)\s*;?'], 'tokens', 'once');
  if ~isempty(tokDef)
    parsed = parseStructArgs(tokDef{1}, localVars, sourceName, '');
    params = mergeParamArrays(params, parsed);
    continue;
  end

  if strcmp(varName, 'default')
    tokInline = regexp(line, '\<readInput\s*\(\s*struct\s*\((.*)\)\s*,', 'tokens', 'once');
    if ~isempty(tokInline)
      parsed = parseStructArgs(tokInline{1}, localVars, sourceName, '');
      params = mergeParamArrays(params, parsed);
    end
  end
end
end


function params = parseStructArgs(argStr, localVars, sourceName, trailingComment)
tokens = splitTopLevelComma(argStr);
params = emptyParamArray();

for k = 1:2:numel(tokens)-1
  rawKey = strtrim(tokens{k});
  rawVal = strtrim(tokens{k+1});

  keyMatch = regexp(rawKey, '^[''"]([^''"]+)[''"]$', 'tokens', 'once');
  if isempty(keyMatch)
    continue;
  end
  keyName = keyMatch{1};

  if isKey(localVars, rawVal)
    rawVal = localVars(rawVal);
  end

  % Skip non-scalar structures and object-handle empty placeholders
  if contains(rawVal, 'struct.empty') || contains(rawVal, 'SimulationParameters(') || ...
      contains(rawVal, 'SimulationParameters.empty') || contains(rawVal, 'Discretizer.empty') || ...
      strcmp(rawVal, 'Grid()') || strcmp(rawVal, 'struct()')
    continue;
  end

  p = classifyValueExpression(keyName, rawVal, sourceName, trailingComment);
  params(end+1) = p; %#ok<AGROW>
end
end


function p = classifyValueExpression(keyName, expr, sourceName, commentText)
p = emptyParam();
p.name = string(keyName);
p.source = string(sourceName);

if ~isempty(commentText)
  p.description = string(commentText);
else
  p.description = sprintf('Parameter of %s', sourceName);
end

exprTrim = strtrim(expr);

if strcmp(exprTrim, 'string.empty')
  p.required = true;
  p.type = "string";
  p.defaultValue = "";
  p.defaultText = "none (compulsory)";
  p.allowedText = "string";
  return;
end

if strcmp(exprTrim, 'double.empty') || strcmp(exprTrim, '[]')
  p.required = true;
  p.type = "double";
  p.defaultValue = double.empty;
  p.defaultText = "none (compulsory)";
  p.allowedText = "real scalar";
  return;
end

if strcmp(exprTrim, 'missing')
  p.required = false;
  p.type = "double";
  p.defaultValue = double.empty;
  p.defaultText = "none (optional)";
  p.allowedText = "real scalar (optional)";
  return;
end

if contains(exprTrim, ':') && (contains(exprTrim, 'nTag') || contains(exprTrim, 'Tag'))
  p.required = false;
  p.type = "vector";
  p.defaultValue = "all";
  p.defaultText = "all region tags";
  p.allowedText = "comma-separated region IDs (e.g. 1,2) or 'all'";
  return;
end

if strcmp(exprTrim, 'true') || strcmp(exprTrim, 'false')
  p.required = false;
  p.type = "flag01";
  p.defaultValue = double(strcmp(exprTrim, 'true'));
  p.defaultText = string(p.defaultValue);
  p.allowedText = "0 or 1 (boolean flag)";
  return;
end

strTok = regexp(exprTrim, '^[''"]([^''"]*)[''"]$', 'tokens', 'once');
if ~isempty(strTok)
  p.required = false;
  p.type = "string";
  p.defaultValue = string(strTok{1});
  p.defaultText = p.defaultValue;
  p.allowedText = "string";
  return;
end

numVal = str2double(exprTrim);
if ~isnan(numVal)
  p.required = false;
  p.type = "double";
  p.defaultValue = numVal;
  p.defaultText = string(exprTrim);
  p.allowedText = "real scalar";
  return;
end

p.required = false;
p.type = "string";
p.defaultValue = string(exprTrim);
p.defaultText = string(exprTrim);
p.allowedText = "value";
end


function branchLines = extractBranchLines(cleanLines, tagName)
% Extract lines inside `if ...<tagName>...` up to matching elseif/end.
branchLines = {};
capturing = false;
depth = 0;
for i = 1:numel(cleanLines)
  line = cleanLines{i};
  if ~capturing
    if contains(line, tagName) && ~isempty(regexp(line, '^\s*(if|elseif)\>', 'once'))
      capturing = true;
      depth = 1;
    end
  else
    if ~isempty(regexp(line, '^\s*(elseif|else)\>', 'once')) && depth == 1
      break;
    end
    sanitized = stripStringsAndBrackets(line);
    words = regexp(sanitized, '\<[A-Za-z_]\w*\>', 'match');
    for w = 1:numel(words)
      if any(strcmp(words{w}, {'if', 'for', 'while', 'switch', 'try'}))
        depth = depth + 1;
      elseif strcmp(words{w}, 'end')
        depth = depth - 1;
      end
    end
    if depth <= 0
      break;
    end
    branchLines{end+1, 1} = line; %#ok<AGROW>
  end
end
end


function branchLines = extractCaseBranchLines(cleanLines, caseVal)
% Extract lines inside `case '<caseVal>'` up to next case/otherwise/end.
branchLines = {};
capturing = false;
depth = 0;
pat = ['^\s*case\s+[''"]' caseVal '[''"]'];
for i = 1:numel(cleanLines)
  line = cleanLines{i};
  if ~capturing
    if ~isempty(regexp(line, pat, 'once'))
      capturing = true;
      depth = 1;
    end
  else
    if depth == 1 && ~isempty(regexp(line, '^\s*(case|otherwise)\>', 'once'))
      break;
    end
    sanitized = stripStringsAndBrackets(line);
    words = regexp(sanitized, '\<[A-Za-z_]\w*\>', 'match');
    for w = 1:numel(words)
      if any(strcmp(words{w}, {'if', 'for', 'while', 'switch', 'try'}))
        depth = depth + 1;
      elseif strcmp(words{w}, 'end')
        depth = depth - 1;
      end
    end
    if depth <= 0
      break;
    end
    branchLines{end+1, 1} = line; %#ok<AGROW>
  end
end
end


function overrides = extractSuperclassOverrides(regBody, methodName)
overrides = containers.Map();
[cleanLines, ~] = collapseContinuations(regBody);

for i = 1:numel(cleanLines)
  tok = regexp(cleanLines{i}, [methodName '@\w+\s*\((.*)\)'], 'tokens', 'once');
  if isempty(tok)
    continue;
  end
  args = splitTopLevelComma(tok{1});
  idx = 1;
  while idx <= numel(args) && isempty(regexp(strtrim(args{idx}), '^[''"]', 'once'))
    idx = idx + 1;
  end
  for k = idx:2:numel(args)-1
    kMatch = regexp(strtrim(args{k}), '^[''"]([^''"]+)[''"]$', 'tokens', 'once');
    vMatch = regexp(strtrim(args{k+1}), '^[''"]([^''"]*)[''"]$', 'tokens', 'once');
    if ~isempty(kMatch) && ~isempty(vMatch)
      overrides(kMatch{1}) = string(vMatch{1});
    end
  end
end
end


function helperNames = extractVararginHelperCalls(regBody, clsName)
helperNames = {};
joined = strjoin(regBody, newline);
mc = meta.class.fromName(clsName);
if isempty(mc)
  return;
end
methodNames = {mc.MethodList.Name};

tok = regexp(joined, '(\w+)\s*\([^)]*varargin\{:\}[^)]*\)', 'tokens');
for i = 1:numel(tok)
  mName = tok{i}{1};
  if ~strcmp(mName, 'registerSolver') && ~strcmp(mName, 'readInput') && ...
      any(strcmp(mName, methodNames))
    helperNames{end+1} = mName; %#ok<AGROW>
  end
end
helperNames = unique(helperNames);
end


function solvers = extractDelegatedSolvers(regBody)
solvers = {};
joined = strjoin(regBody, newline);
tok = regexp(joined, '=\s*(\w+)\s*\(\s*obj\.domain\s*\)', 'tokens');
for i = 1:numel(tok)
  cand = tok{i}{1};
  if exist(cand, 'class') == 8
    solvers{end+1} = cand; %#ok<AGROW>
  end
end
solvers = unique(solvers);
end


function params = applyOverrides(params, overrideMap)
if isempty(overrideMap) || overrideMap.Count == 0
  return;
end
for i = 1:numel(params)
  k = char(params(i).name);
  if isKey(overrideMap, k)
    val = overrideMap(k);
    params(i).defaultValue = val;
    params(i).defaultText = string(val);
  end
end
end


function out = mergeParamArrays(baseArr, newArr)
out = baseArr;
for i = 1:numel(newArr)
  if isempty(out)
    out(end+1) = newArr(i); %#ok<AGROW>
  else
    idx = find([out.name] == newArr(i).name, 1);
    if isempty(idx)
      out(end+1) = newArr(i); %#ok<AGROW>
    else
      out(idx) = newArr(i);
    end
  end
end
end


function p = emptyParam()
p = struct('name', "", 'type', "double", 'defaultValue', [], ...
  'defaultText', "", 'required', false, 'allowedText', "", ...
  'choices', string.empty, 'description', "", 'source', "", 'helpText', "");
end


function arr = emptyParamArray()
arr = repmat(emptyParam(), 1, 0);
end


function lines = readLines(filePath)
txt = fileread(filePath);
lines = regexp(txt, '\r\n|\n|\r', 'split')';
end


function [cleanLines, comments] = collapseContinuations(lines)
cleanLines = {};
comments = {};
buf = '';
commentBuf = '';

for i = 1:numel(lines)
  raw = lines{i};
  cIdx = findPercentOutsideStrings(raw);
  if ~isempty(cIdx)
    codePart = raw(1:cIdx-1);
    cText = strtrim(raw(cIdx+1:end));
    if ~isempty(cText) && isempty(commentBuf)
      commentBuf = cText;
    end
  else
    codePart = raw;
  end

  trimmed = strtrim(codePart);
  if endsWith(trimmed, '...')
    buf = [buf ' ' strtrim(trimmed(1:end-3))]; %#ok<AGROW>
  else
    fullStmt = strtrim([buf ' ' trimmed]);
    if ~isempty(fullStmt)
      cleanLines{end+1} = fullStmt; %#ok<AGROW>
      comments{end+1} = commentBuf; %#ok<AGROW>
    end
    buf = '';
    commentBuf = '';
  end
end
end


function idx = findPercentOutsideStrings(line)
idx = [];
inSingle = false;
inDouble = false;
for k = 1:length(line)
  ch = line(k);
  if ch == '''' && ~inDouble
    inSingle = ~inSingle;
  elseif ch == '"' && ~inSingle
    inDouble = ~inDouble;
  elseif ch == '%' && ~inSingle && ~inDouble
    idx = k;
    return;
  end
end
end


function tokens = splitTopLevelComma(str)
tokens = {};
depth = 0;
inSingle = false;
inDouble = false;
startIdx = 1;

for k = 1:length(str)
  ch = str(k);
  if ch == '''' && ~inDouble
    inSingle = ~inSingle;
  elseif ch == '"' && ~inSingle
    inDouble = ~inDouble;
  elseif ~inSingle && ~inDouble
    if ch == '(' || ch == '[' || ch == '{'
      depth = depth + 1;
    elseif ch == ')' || ch == ']' || ch == '}'
      depth = depth - 1;
    elseif ch == ',' && depth == 0
      tokens{end+1} = str(startIdx:k-1); %#ok<AGROW>
      startIdx = k + 1;
    end
  end
end
tokens{end+1} = str(startIdx:end);
end


function body = extractMethodBody(lines, methodName)
body = {};
startLine = [];
pat = ['^\s*function\>\s*([^=]*=\s*)?\<' methodName '\>'];
for i = 1:numel(lines)
  if ~isempty(regexp(stripComment(lines{i}), pat, 'once'))
    startLine = i;
    break;
  end
end
if isempty(startLine)
  return;
end

depth = 1;
body{end+1,1} = lines{startLine};
openKeywords = {'function', 'if', 'for', 'while', 'switch', 'try', 'parfor', 'spmd'};

for i = startLine+1:numel(lines)
  code = strtrim(stripComment(lines{i}));
  body{end+1,1} = lines{i}; %#ok<AGROW>
  if isempty(code)
    continue;
  end

  sanitized = stripStringsAndBrackets(code);
  words = regexp(sanitized, '\<[A-Za-z_]\w*\>', 'match');
  for w = 1:numel(words)
    kw = words{w};
    if any(strcmp(kw, openKeywords))
      depth = depth + 1;
    elseif strcmp(kw, 'end')
      depth = depth - 1;
    end
  end
  if depth <= 0
    break;
  end
end
end


function out = stripStringsAndBrackets(line)
out = blanks(length(line));
inSingle = false;
inDouble = false;
parenDepth = 0;

for k = 1:length(line)
  ch = line(k);
  if ch == '''' && ~inDouble
    inSingle = ~inSingle;
  elseif ch == '"' && ~inSingle
    inDouble = ~inDouble;
  elseif ~inSingle && ~inDouble
    if ch == '(' || ch == '[' || ch == '{'
      parenDepth = parenDepth + 1;
    elseif ch == ')' || ch == ']' || ch == '}'
      parenDepth = max(0, parenDepth - 1);
    elseif parenDepth == 0
      out(k) = ch;
    end
  end
end
end


function val = findPropertyDefault(lines, propName)
val = [];
for i = 1:numel(lines)
  code = strtrim(stripComment(lines{i}));
  tok = regexp(code, ['^' propName '\s*=\s*([^;%]+)'], 'tokens', 'once');
  if ~isempty(tok)
    expr = strtrim(tok{1});
    if strcmp(expr, 'true')
      val = true;
      return;
    elseif strcmp(expr, 'false')
      val = false;
      return;
    end
    sTok = regexp(expr, '^[''"]([^''"]*)[''"]$', 'tokens', 'once');
    if ~isempty(sTok)
      val = string(sTok{1});
      return;
    end
    nVal = str2double(expr);
    if ~isnan(nVal)
      val = nVal;
      return;
    end
  end
end
end


function desc = extractLeadingComment(lines, lineIdx)
desc = '';
for k = lineIdx-1:-1:max(1, lineIdx-3)
  t = strtrim(lines{k});
  if startsWith(t, '%')
    desc = strtrim(regexprep(t, '^%+\s*', ''));
    return;
  elseif ~isempty(t)
    return;
  end
end
end


function code = stripComment(line)
idx = findPercentOutsideStrings(line);
if isempty(idx)
  code = line;
else
  code = line(1:idx-1);
end
end


function chain = getInheritanceChain(mc)
chain = {};
if ~isempty(mc.SuperclassList)
  for i = 1:numel(mc.SuperclassList)
    subChain = getInheritanceChain(mc.SuperclassList(i));
    chain = [chain, subChain]; %#ok<AGROW>
  end
end
chain{end+1} = mc.Name;
chain = unique(chain, 'stable');
end


function yes = inheritsFrom(mc, superName)
yes = false;
if isempty(mc.SuperclassList)
  return;
end
for i = 1:numel(mc.SuperclassList)
  s = mc.SuperclassList(i);
  if strcmp(s.Name, superName) || inheritsFrom(s, superName)
    yes = true;
    return;
  end
end
end


function yes = hasStaticMethod(mc, mName)
idx = strcmp({mc.MethodList.Name}, mName);
yes = any(idx) && mc.MethodList(find(idx, 1)).Static;
end


function files = findClassFiles(rootFolder)
files = findMFilesMatching(rootFolder, true);
end


function files = findMFiles(rootFolder)
files = findMFilesMatching(rootFolder, false);
end


function files = findMFilesMatching(rootFolder, onlyClasses)
files = {};
items = dir(rootFolder);
for i = 1:numel(items)
  name = items(i).name;
  if name(1) == '.'
    continue;
  end
  fullP = fullfile(rootFolder, name);
  if items(i).isdir
    files = [files, findMFilesMatching(fullP, onlyClasses)]; %#ok<AGROW>
  elseif endsWith(name, '.m')
    if ~onlyClasses || contains(fileread(fullP), 'classdef')
      files{end+1} = fullP; %#ok<AGROW>
    end
  end
end
end
