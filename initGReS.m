function initGReS(outputFlag)
% List of folders to check
foldList = ["Code",...
            "Utilities",...
            "ThidPartyLibs/Chronos_Lab"];

if nargin == 0
  outputFlag = true;
end

% Check if any of the folders is already on the MATLAB path
needRestore = false;
for k = 1:numel(foldList)
  if contains(path, fullfile('GReS',foldList(k)))
    needRestore = true;
    break;
  end
end

% Restore path only if needed
if needRestore
  if outputFlag
    fprintf('Restoring default paths \n');
  end
  restoredefaultpath;
end

% set the root to the gres directory
gres_root = fileparts(mfilename('fullpath'));
setappdata(0,'gres_root', gres_root);

% Add GReS paths
for k = 1:numel(foldList)
  addpath(genpath(fullfile(gres_root, foldList(k))));
end

setappdata(0,'gresLog', Logger());
if outputFlag
  gresLog().welcomeMsg()
end


% initialize MRST if available
mrstPath = fullfile(gres_root,'ThirdPartyLibs','MRST');
if isfile(fullfile(mrstPath,'startup.m'))
  cd(mrstPath)
  [~] = evalc('startup');
else
  warning(['MRST not active. To get MRST, run:' ...
    ' git submodule update --init --recursive ThirdPartyLibs/MRST in the GReS root directory'])
end

cd(gres_root)


end