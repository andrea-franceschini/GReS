function compileAll(varargin)

if isempty(varargin)
  compileList = ["gres","chronos"];
else
  compileList = string(varargin);
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
chronosDir = fullfile(gres_root,'ThirdPartyLibs','Chronos_Lab');
isChronosReady = isfile(fullfile(chronosDir,'compileAll.m'));

if isChronosReady
   cd(chronosDir);
   compileAll
   cd(home)
end

end



end
