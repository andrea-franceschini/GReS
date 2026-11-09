function compileAll(varargin)
%COMPILEALL Build GReS MEX files and optional installed dependencies.
%   compileAll('gmsh') builds only the Gmsh MEX and requires setupGmsh first.

if isempty(varargin)
  compileList = ["gres","chronos","gmsh"];
  requireGmsh = false;
else
  compileList = string(varargin);
  requireGmsh = ismember("gmsh",compileList);
end

clc;
close all;

if ismember("gres",compileList)

list = {'Code/grid',...
        'Code/output', ...
        'Code/elements',...
        'Utilities/ComputationalGeometry'};

home = pwd;
for folder = list
  fprintf('Compiling MEX files in %s\n', folder{1});
  cd(folder{1});
  compile
  fprintf('\n\n')
  cd(home);
end

end

if ismember("chronos",compileList)

% compiling Chronos if available
chronosDir = fullfile(gres_root,'ThirdPartyLibs','ChronosLab');
isChronosReady = isfile(fullfile(chronosDir,'compileAll.m'));

if isChronosReady
   cd(chronosDir);
   compileAll
   cd(gres_root)
end

end

if ismember("gmsh",compileList)
  [sdk, ready] = gmshSdk();
  if ready
    mexDir = fullfile(gres_root,'Code','grid','mex');
    args = ["-O", "-outdir", string(mexDir), ...
      "-I" + sdk.includeDir, fullfile(mexDir,'gmshMex.cpp'), sdk.linkFile];
    if ~ispc
      % Let the loader find the SDK without changing library search variables.
      args(end+1) = 'LDFLAGS=$LDFLAGS -Wl,-rpath,"' + sdk.libDir + '"';
    end
    fprintf('Compiling Gmsh MEX with %s\n', sdk.name);
    mex(args{:});
    addGmshRuntimePath(sdk);
  elseif requireGmsh
    error('GReS:GmshSDKUnavailable', ...
      'Gmsh SDK unavailable for %s. Run setupGmsh to enable Gmsh support.', computer);
  else
    fprintf('Gmsh SDK unavailable for %s. Skipping Gmsh interface.\n', computer);
    fprintf('Run setupGmsh to enable Gmsh support.\n');
  end
end

end
