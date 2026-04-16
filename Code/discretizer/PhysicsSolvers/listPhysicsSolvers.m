function listPhysicsSolvers()

rootFolder = fileparts(mfilename('fullpath'));
classFiles = findClassesRecursive(rootFolder);

fprintf('=====================================================\n');
fprintf('   List of Physics Solvers and registered variables  \n');
fprintf('=====================================================\n');

for i = 1:numel(classFiles)
  className = classFiles{i};

  % Load meta class
  mc = meta.class.fromName(className);
  if isempty(mc)
    continue;
  end

  % Must inherit from PhysicsSolver
  if ~inheritsFrom(mc, 'PhysicsSolver')
    continue;
  end

  % Must implement a static getField()
  idx = strcmp({mc.MethodList.Name}, 'getField');
  if ~any(idx)
    continue;
  end

  m = mc.MethodList(idx);
  if ~m.Static
    continue;
  end

  % Call getField
  try
    fields = feval(className + ".getField");
  catch
    continue;
  end

  % Print formatted output
  fprintf(' \n Solver Class: %s\n', className);
  fprintf('   Registers fields:\n');

  if isempty(fields)
    fprintf('      (none)\n\n');
  elseif iscell(fields)
    for k = 1:numel(fields)
      fprintf('      - %s\n', fields{k});
    end
    fprintf('\n');
  else
    fprintf('      - %s\n', string(fields));
  end
end

fprintf('===============================================\n');
end


function classList = findClassesRecursive(rootFolder)

classList = {};
items = dir(rootFolder);

for i = 1:numel(items)
  name = items(i).name;

  % skip .  and ..
  if name(1) == '.'
    continue;
  end

  pathItem = fullfile(rootFolder, name);

  if items(i).isdir
    sub = findClassesRecursive(pathItem);
    classList = [classList, sub]; %#ok<AGROW>
    continue;
  end

  if endsWith(name, '.m')
    fileText = fileread(pathItem);
    if contains(fileText, 'classdef')
      [~, className] = fileparts(name);
      classList{end+1} = className; %#ok<AGROW>
    end
  end
end
end



function yes = inheritsFrom(mc, superName)
yes = false;
for sc = mc.SuperclassList
  if strcmp(sc.Name, superName)
    yes = true;
    return;
  else
    % Recursively check parent classes
    yes = inheritsFrom(sc, superName);
    if yes
      return;
    end
  end
end
end