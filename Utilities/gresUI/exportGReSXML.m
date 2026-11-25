function writtenFiles = exportGReSXML(config, outputPath, mode)
%EXPORTGRESXML Serialize gresUI selections to GReS-compatible XML file(s).
%
%   writtenFiles = exportGReSXML(config, outputPath, mode)
%
%   INPUTS:
%     config     : struct with .simParams, .linearSolver, .domains (array of
%                  domain structs with .solver and optional .geometry,
%                  .materials, .boundaryConditions), plus optional
%                  .interfaces and .output
%     outputPath : file path (unified mode) or directory/base path (modular mode)
%     mode       : "unified" or "modular"
%
%   OUTPUT:
%     writtenFiles : string array of generated XML file path(s)

mode = lower(string(mode));
validateattributes(mode, {'string', 'char'}, {'scalartext'});

simXml = buildSimParamsXml(config.simParams, config.linearSolver, 0);
domainsCfg = config.domains;
nDom = numel(domainsCfg);
interfacesCfg = getFieldSafe(config, 'interfaces');
hasInterf = ~isempty(interfacesCfg) && isstruct(interfacesCfg);
interfXml = buildInterfacesXml(interfacesCfg, 0);
outXml = buildOutputXml(getFieldSafe(config, 'output'), 0);

switch mode
  case "unified"
    blocks = string.empty;
    blocks(end+1) = indentBlock(simXml, 2);

    for d = 1:nDom
      dom = domainsCfg(d);
      solverXml = buildSolverXml(dom.solver, 0);
      geomXml = buildGeometryXml(getFieldSafe(dom, 'geometry'), 0);
      matXml = buildMaterialsXml(getFieldSafe(dom, 'materials'), 0);
      bcXml = buildBoundariesXml(getFieldSafe(dom, 'boundaryConditions'), 0);
      hasDomainChildren = strlength(geomXml) > 0 || strlength(matXml) > 0 || strlength(bcXml) > 0;

      if nDom == 1 && ~hasInterf && ~hasDomainChildren
        blocks(end+1) = indentBlock(solverXml, 2); %#ok<AGROW>
      else
        domChildren = string.empty;
        if strlength(geomXml) > 0, domChildren(end+1) = indentBlock(geomXml, 4); end %#ok<AGROW>
        if strlength(matXml) > 0, domChildren(end+1) = indentBlock(matXml, 4); end %#ok<AGROW>
        if strlength(bcXml) > 0, domChildren(end+1) = indentBlock(bcXml, 4); end %#ok<AGROW>
        domChildren(end+1) = indentBlock(solverXml, 4); %#ok<AGROW>
        blocks(end+1) = "  <Domain>" + newline + strjoin(domChildren, [newline newline]) + newline + "  </Domain>"; %#ok<AGROW>
      end
    end

    if strlength(interfXml) > 0
      blocks(end+1) = indentBlock(interfXml, 2);
    end

    if strlength(outXml) > 0
      blocks(end+1) = indentBlock(outXml, 2);
    end

    xmlText = "<?xml version=""1.0"" encoding=""UTF-8""?>" + newline + ...
      "<Problem>" + newline + strjoin(blocks, [newline newline]) + newline + "</Problem>" + newline;
    writeTextFile(outputPath, xmlText);
    writtenFiles = string(outputPath);

  case "modular"
    if isfolder(outputPath)
      outDir = char(outputPath);
    else
      [outDir, ~, ~] = fileparts(char(outputPath));
      if isempty(outDir)
        outDir = pwd;
      end
      if ~isfolder(outDir)
        mkdir(outDir);
      end
    end

    header = "<?xml version=""1.0"" encoding=""UTF-8""?>" + newline;
    writtenFiles = string.empty;

    simPath = fullfile(outDir, 'simparams.xml');
    writeTextFile(simPath, header + simXml + newline);
    writtenFiles(end+1) = string(simPath);

    probBlocks = string.empty;
    probBlocks(end+1) = "  <SimulationParameters fileName=""simparams.xml""/>";

    for d = 1:nDom
      dom = domainsCfg(d);
      solverXml = buildSolverXml(dom.solver, 0);
      geomXml = buildGeometryXml(getFieldSafe(dom, 'geometry'), 0);
      matXml = buildMaterialsXml(getFieldSafe(dom, 'materials'), 0);
      bcXml = buildBoundariesXml(getFieldSafe(dom, 'boundaryConditions'), 0);

      if nDom == 1
        prefix = '';
        domIndent = 4;
      else
        prefix = sprintf('domain%d_', d);
        domIndent = 2;
      end

      solName = [prefix 'solver.xml'];
      solPath = fullfile(outDir, solName);
      writeTextFile(solPath, header + solverXml + newline);
      writtenFiles(end+1) = string(solPath); %#ok<AGROW>

      domRefs = string.empty;
      if strlength(geomXml) > 0
        domRefs(end+1) = indentBlock(geomXml, domIndent); %#ok<AGROW>
      end

      if strlength(matXml) > 0
        matName = [prefix 'materials.xml'];
        matPath = fullfile(outDir, matName);
        writeTextFile(matPath, header + matXml + newline);
        writtenFiles(end+1) = string(matPath); %#ok<AGROW>
        domRefs(end+1) = repmat(' ', 1, domIndent) + "<Materials fileName=""" + string(matName) + """/>"; %#ok<AGROW>
      end

      if strlength(bcXml) > 0
        bcName = [prefix 'boundaries.xml'];
        bcPath = fullfile(outDir, bcName);
        writeTextFile(bcPath, header + bcXml + newline);
        writtenFiles(end+1) = string(bcPath); %#ok<AGROW>
        domRefs(end+1) = repmat(' ', 1, domIndent) + "<BoundaryConditions fileName=""" + string(bcName) + """/>"; %#ok<AGROW>
      end

      domRefs(end+1) = repmat(' ', 1, domIndent) + "<Solver fileName=""" + string(solName) + """/>"; %#ok<AGROW>

      if nDom == 1
        probBlocks(end+1) = "  <Domain>" + newline + strjoin(domRefs, newline) + newline + "  </Domain>"; %#ok<AGROW>
      else
        domFileName = sprintf('domain%d.xml', d);
        domPath = fullfile(outDir, domFileName);
        domContent = header + "<Domain>" + newline + strjoin(domRefs, newline) + newline + "</Domain>" + newline;
        writeTextFile(domPath, domContent);
        writtenFiles(end+1) = string(domPath); %#ok<AGROW>
        probBlocks(end+1) = "  <Domain fileName=""" + string(domFileName) + """/>"; %#ok<AGROW>
      end
    end

    if strlength(interfXml) > 0
      interfPath = fullfile(outDir, 'interfaces.xml');
      writeTextFile(interfPath, header + interfXml + newline);
      writtenFiles(end+1) = string(interfPath);
      probBlocks(end+1) = "  <Interface fileName=""interfaces.xml""/>";
    end

    if strlength(outXml) > 0
      outPath = fullfile(outDir, 'output.xml');
      writeTextFile(outPath, header + outXml + newline);
      writtenFiles(end+1) = string(outPath);
      probBlocks(end+1) = "  <Output fileName=""output.xml""/>";
    end

    probPath = fullfile(outDir, 'problem.xml');
    probText = header + "<Problem>" + newline + strjoin(probBlocks, newline) + newline + "</Problem>" + newline;
    writeTextFile(probPath, probText);
    writtenFiles(end+1) = string(probPath);

  otherwise
    error('exportGReSXML:InvalidMode', 'Unsupported export mode "%s". Expected "unified" or "modular".', mode);
end
end


function val = getFieldSafe(s, fName)
if isstruct(s) && isfield(s, fName)
  val = s.(fName);
else
  val = [];
end
end


function xml = buildSimParamsXml(simParams, linSolver, indent)
pad = repmat(' ', 1, indent);
attrs = formatAttributes(simParams, indent + 2);
linXml = buildLinearSolverXml(linSolver, indent + 2);

if strlength(linXml) == 0
  if strlength(attrs) == 0
    xml = pad + "<SimulationParameters/>";
  else
    xml = pad + "<SimulationParameters" + newline + attrs + "/>";
  end
else
  if strlength(attrs) == 0
    xml = pad + "<SimulationParameters>" + newline + linXml + newline + pad + "</SimulationParameters>";
  else
    xml = pad + "<SimulationParameters" + newline + attrs + ">" + newline + ...
      linXml + newline + pad + "</SimulationParameters>";
  end
end
end


function xml = buildLinearSolverXml(linSolver, indent)
if isempty(linSolver) || ~isstruct(linSolver)
  xml = "";
  return;
end

pad = repmat(' ', 1, indent);
topAttrs = "";
if isfield(linSolver, 'topLevel') && isstruct(linSolver.topLevel) && ~isempty(fieldnames(linSolver.topLevel))
  topAttrs = formatAttributes(linSolver.topLevel, indent + 2);
end

children = string.empty;
physKeys = ["Flow", "Mechanics"];
for k = 1:numel(physKeys)
  pk = physKeys(k);
  if isfield(linSolver, pk) && isstruct(linSolver.(pk))
    subXml = buildChronosPhysicsXml(pk, linSolver.(pk), indent + 2);
    if strlength(subXml) > 0
      children(end+1) = subXml; %#ok<AGROW>
    end
  end
end

if strlength(topAttrs) == 0 && isempty(children)
  xml = "";
  return;
end

if isempty(children)
  xml = pad + "<LinearSolver" + newline + topAttrs + "/>";
else
  openTag = pad + "<LinearSolver";
  if strlength(topAttrs) > 0
    openTag = openTag + newline + topAttrs;
  end
  xml = openTag + ">" + newline + strjoin(children, newline) + newline + pad + "</LinearSolver>";
end
end


function xml = buildChronosPhysicsXml(physName, physStruct, indent)
pad = repmat(' ', 1, indent);
secNames = fieldnames(physStruct);
secNodes = string.empty;

for i = 1:numel(secNames)
  sec = secNames{i};
  secData = physStruct.(sec);
  if ~isstruct(secData) || isempty(fieldnames(secData))
    continue;
  end
  attrs = formatAttributes(secData, indent + 4);
  if strlength(attrs) > 0
    secNodes(end+1) = repmat(' ', 1, indent + 2) + "<" + string(sec) + newline + attrs + "/>"; %#ok<AGROW>
  end
end

if isempty(secNodes)
  xml = "";
else
  xml = pad + "<" + physName + ">" + newline + ...
    strjoin(secNodes, newline) + newline + ...
    pad + "</" + physName + ">";
end
end


function xml = buildGeometryXml(geomCfg, indent)
if isempty(geomCfg) || ~isstruct(geomCfg) || ~isfield(geomCfg, 'type')
  xml = "";
  return;
end
pad = repmat(' ', 1, indent);
innerPad = repmat(' ', 1, indent + 2);
gType = string(geomCfg.type);
if gType == "meshFile"
  xml = pad + "<Geometry meshFile=""" + escapeXml(string(geomCfg.params.meshFile)) + """/>";
else
  attrs = formatAttributes(geomCfg.params, indent + 4);
  xml = pad + "<Geometry>" + newline + ...
    innerPad + "<" + gType + newline + attrs + "/>" + newline + ...
    pad + "</Geometry>";
end
end


function xml = buildMaterialsXml(matCfg, indent)
if isempty(matCfg) || ~isstruct(matCfg)
  xml = "";
  return;
end
pad = repmat(' ', 1, indent);
innerPad = repmat(' ', 1, indent + 2);
children = string.empty;

if isfield(matCfg, 'Fluid') && isstruct(matCfg.Fluid) && ~isempty(fieldnames(matCfg.Fluid))
  fAttrs = formatAttributes(matCfg.Fluid, indent + 4);
  if strlength(fAttrs) > 0
    children(end+1) = innerPad + "<Fluid" + newline + fAttrs + "/>";
  else
    children(end+1) = innerPad + "<Fluid/>";
  end
end

if isfield(matCfg, 'Solids') && ~isempty(matCfg.Solids)
  for i = 1:numel(matCfg.Solids)
    s = matCfg.Solids(i);
    sAttrs = formatAttributes(s.params, indent + 4);
    subNodes = string.empty;

    if isfield(s, 'constitutive') && isstruct(s.constitutive) && isfield(s.constitutive, 'name') && strlength(string(s.constitutive.name)) > 0
      cName = string(s.constitutive.name);
      cAttrs = formatAttributes(s.constitutive.params, indent + 8);
      if strlength(cAttrs) > 0
        lawNode = repmat(' ', 1, indent + 6) + "<" + cName + newline + cAttrs + "/>";
      else
        lawNode = repmat(' ', 1, indent + 6) + "<" + cName + "/>";
      end
      subNodes(end+1) = repmat(' ', 1, indent + 4) + "<Constitutive>" + newline + ...
        lawNode + newline + repmat(' ', 1, indent + 4) + "</Constitutive>"; %#ok<AGROW>
    end

    if isfield(s, 'porousRock') && isstruct(s.porousRock)
      prAttrs = formatAttributes(s.porousRock.params, indent + 6);
      hasCurves = isfield(s.porousRock, 'curves') && isstruct(s.porousRock.curves) && ~isempty(fieldnames(s.porousRock.curves));
      if hasCurves
        cvAttrs = formatAttributes(s.porousRock.curves, indent + 8);
        cvNode = repmat(' ', 1, indent + 6) + "<Curves" + newline + cvAttrs + "/>";
        if strlength(prAttrs) > 0
          subNodes(end+1) = repmat(' ', 1, indent + 4) + "<PorousRock" + newline + prAttrs + ">" + newline + ...
            cvNode + newline + repmat(' ', 1, indent + 4) + "</PorousRock>"; %#ok<AGROW>
        else
          subNodes(end+1) = repmat(' ', 1, indent + 4) + "<PorousRock>" + newline + ...
            cvNode + newline + repmat(' ', 1, indent + 4) + "</PorousRock>"; %#ok<AGROW>
        end
      else
        if strlength(prAttrs) > 0
          subNodes(end+1) = repmat(' ', 1, indent + 4) + "<PorousRock" + newline + prAttrs + "/>"; %#ok<AGROW>
        else
          subNodes(end+1) = repmat(' ', 1, indent + 4) + "<PorousRock/>"; %#ok<AGROW>
        end
      end
    end

    if isempty(subNodes)
      children(end+1) = innerPad + "<Solid" + newline + sAttrs + "/>"; %#ok<AGROW>
    else
      children(end+1) = innerPad + "<Solid" + newline + sAttrs + ">" + newline + ...
        strjoin(subNodes, newline) + newline + innerPad + "</Solid>"; %#ok<AGROW>
    end
  end
end

if isempty(children)
  xml = "";
else
  xml = pad + "<Materials>" + newline + strjoin(children, newline) + newline + pad + "</Materials>";
end
end


function xml = buildBoundariesXml(bcs, indent)
if isempty(bcs) || ~isstruct(bcs)
  xml = "";
  return;
end
pad = repmat(' ', 1, indent);
innerPad = repmat(' ', 1, indent + 2);
bcNodes = string.empty;

for i = 1:numel(bcs)
  bc = bcs(i);
  attrs = formatAttributes(bc.params, indent + 4);
  evNodes = string.empty;
  if isfield(bc, 'events') && ~isempty(bc.events)
    for e = 1:numel(bc.events)
      evAttrs = formatAttributes(bc.events(e), indent + 6);
      evNodes(end+1) = repmat(' ', 1, indent + 4) + "<BCevent" + newline + evAttrs + "/>"; %#ok<AGROW>
    end
  end
  if isempty(evNodes)
    bcNodes(end+1) = innerPad + "<BC" + newline + attrs + "/>"; %#ok<AGROW>
  else
    bcNodes(end+1) = innerPad + "<BC" + newline + attrs + ">" + newline + ...
      strjoin(evNodes, newline) + newline + innerPad + "</BC>"; %#ok<AGROW>
  end
end

if isempty(bcNodes)
  xml = "";
else
  xml = pad + "<BoundaryConditions>" + newline + strjoin(bcNodes, newline) + newline + pad + "</BoundaryConditions>";
end
end


function xml = buildInterfacesXml(ifaces, indent)
if isempty(ifaces) || ~isstruct(ifaces)
  xml = "";
  return;
end
pad = repmat(' ', 1, indent);
innerPad = repmat(' ', 1, indent + 2);
modelNodes = string.empty;

for i = 1:numel(ifaces)
  ifCfg = ifaces(i);
  mName = string(ifCfg.type);
  attrs = formatAttributes(getFieldSafe(ifCfg, 'params'), indent + 4);
  subNodes = string.empty;

  quad = getFieldSafe(ifCfg, 'quadrature');
  if isstruct(quad) && isfield(quad, 'type')
    qStruct = struct('type', quad.type);
    if isfield(quad, 'params') && isstruct(quad.params)
      qFn = fieldnames(quad.params);
      for qk = 1:numel(qFn)
        qStruct.(qFn{qk}) = quad.params.(qFn{qk});
      end
    end
    qAttrs = formatAttributes(qStruct, indent + 6);
    subNodes(end+1) = repmat(' ', 1, indent + 4) + "<Quadrature" + newline + qAttrs + "/>"; %#ok<AGROW>
  end

  coulomb = getFieldSafe(ifCfg, 'coulomb');
  if isstruct(coulomb) && ~isempty(fieldnames(coulomb))
    cAttrs = formatAttributes(coulomb, indent + 6);
    subNodes(end+1) = repmat(' ', 1, indent + 4) + "<Coulomb" + newline + cAttrs + "/>"; %#ok<AGROW>
  end

  as = getFieldSafe(ifCfg, 'activeSet');
  if isstruct(as)
    asAttrs = formatAttributes(getFieldSafe(as, 'params'), indent + 6);
    tolStruct = getFieldSafe(as, 'tolerances');
    hasTol = isstruct(tolStruct) && ~isempty(fieldnames(tolStruct));
    if hasTol
      tAttrs = formatAttributes(tolStruct, indent + 8);
      tNode = repmat(' ', 1, indent + 6) + "<Tolerances" + newline + tAttrs + "/>";
      openAs = repmat(' ', 1, indent + 4) + "<ActiveSet";
      if strlength(asAttrs) > 0
        openAs = openAs + newline + asAttrs;
      end
      subNodes(end+1) = openAs + ">" + newline + tNode + newline + ...
        repmat(' ', 1, indent + 4) + "</ActiveSet>"; %#ok<AGROW>
    elseif strlength(asAttrs) > 0
      subNodes(end+1) = repmat(' ', 1, indent + 4) + "<ActiveSet" + newline + asAttrs + "/>"; %#ok<AGROW>
    else
      subNodes(end+1) = repmat(' ', 1, indent + 4) + "<ActiveSet/>"; %#ok<AGROW>
    end
  end

  if isempty(subNodes)
    if strlength(attrs) == 0
      modelNodes(end+1) = innerPad + "<" + mName + "/>"; %#ok<AGROW>
    else
      modelNodes(end+1) = innerPad + "<" + mName + newline + attrs + "/>"; %#ok<AGROW>
    end
  else
    openM = innerPad + "<" + mName;
    if strlength(attrs) > 0
      openM = openM + newline + attrs;
    end
    modelNodes(end+1) = openM + ">" + newline + strjoin(subNodes, newline) + newline + ...
      innerPad + "</" + mName + ">"; %#ok<AGROW>
  end
end

if isempty(modelNodes)
  xml = "";
else
  xml = pad + "<Interface>" + newline + strjoin(modelNodes, newline) + newline + pad + "</Interface>";
end
end


function xml = buildOutputXml(outCfg, indent)
if isempty(outCfg) || ~isstruct(outCfg)
  xml = "";
  return;
end
pad = repmat(' ', 1, indent);
attrs = "";
if isfield(outCfg, 'params') && isstruct(outCfg.params)
  attrs = formatAttributes(outCfg.params, indent + 2);
end
hasRestart = isfield(outCfg, 'restart') && isstruct(outCfg.restart) && ~isempty(fieldnames(outCfg.restart));

if strlength(attrs) == 0 && ~hasRestart
  xml = "";
  return;
end

if ~hasRestart
  xml = pad + "<Output" + newline + attrs + "/>";
else
  rAttrs = formatAttributes(outCfg.restart, indent + 4);
  rNode = repmat(' ', 1, indent + 2) + "<restart" + newline + rAttrs + "/>";
  if strlength(attrs) == 0
    xml = pad + "<Output>" + newline + rNode + newline + pad + "</Output>";
  else
    xml = pad + "<Output" + newline + attrs + ">" + newline + rNode + newline + pad + "</Output>";
  end
end
end


function xml = buildSolverXml(solverCfg, indent)
pad = repmat(' ', 1, indent);
innerPad = repmat(' ', 1, indent + 2);

sName = string(solverCfg.name);
attrs = "";
if isfield(solverCfg, 'params') && isstruct(solverCfg.params) && ~isempty(fieldnames(solverCfg.params))
  attrs = formatAttributes(solverCfg.params, indent + 4);
end

children = string.empty;
childKeys = ["Mechanics", "Flow"];
for k = 1:numel(childKeys)
  ck = childKeys(k);
  if isfield(solverCfg, ck) && isstruct(solverCfg.(ck)) && ~isempty(fieldnames(solverCfg.(ck)))
    cAttrs = formatAttributes(solverCfg.(ck), indent + 6);
    if strlength(cAttrs) > 0
      children(end+1) = repmat(' ', 1, indent + 4) + "<" + ck + newline + cAttrs + "/>"; %#ok<AGROW>
    end
  end
end

if isempty(children)
  if strlength(attrs) == 0
    solverNode = innerPad + "<" + sName + "/>";
  else
    solverNode = innerPad + "<" + sName + newline + attrs + "/>";
  end
else
  openTag = innerPad + "<" + sName;
  if strlength(attrs) > 0
    openTag = openTag + newline + attrs;
  end
  solverNode = openTag + ">" + newline + strjoin(children, newline) + newline + innerPad + "</" + sName + ">";
end

xml = pad + "<Solver>" + newline + solverNode + newline + pad + "</Solver>";
end


function attrText = formatAttributes(s, indent)
if isempty(s) || ~isstruct(s)
  attrText = "";
  return;
end

fNames = fieldnames(s);
if isempty(fNames)
  attrText = "";
  return;
end

pad = repmat(' ', 1, indent);
lines = string.empty;
for i = 1:numel(fNames)
  fn = fNames{i};
  valStr = formatXmlValue(s.(fn));
  lines(end+1) = pad + string(fn) + " = """ + escapeXml(valStr) + """"; %#ok<AGROW>
end
attrText = strjoin(lines, newline);
end


function str = formatXmlValue(val)
if islogical(val)
  if isscalar(val)
    str = string(double(val));
  else
    str = strjoin(string(double(val(:))), ",");
  end
elseif isnumeric(val)
  if isscalar(val)
    str = string(sprintf('%.15g', val));
  else
    parts = arrayfun(@(x) sprintf('%.15g', x), val(:)', 'UniformOutput', false);
    str = string(strjoin(parts, ','));
  end
elseif isstring(val) || ischar(val)
  v = string(val);
  if isscalar(v)
    str = v;
  else
    str = strjoin(v(:), ",");
  end
else
  error('exportGReSXML:UnsupportedType', 'Cannot serialize value of class %s.', class(val));
end
end


function out = escapeXml(str)
out = strrep(str, "&", "&amp;");
out = strrep(out, """", "&quot;");
out = strrep(out, "<", "&lt;");
out = strrep(out, ">", "&gt;");
end


function out = indentBlock(block, spaces)
pad = repmat(' ', 1, spaces);
lines = splitlines(block);
out = strjoin(pad + lines, newline);
end


function writeTextFile(filePath, content)
fid = fopen(filePath, 'w', 'n', 'UTF-8');
if fid == -1
  error('exportGReSXML:FileOpenFailed', 'Could not open "%s" for writing.', filePath);
end
cleanupObj = onCleanup(@() fclose(fid));
fwrite(fid, char(content), 'char');
end
