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
  compileGmsh(requireGmsh);
end

end
