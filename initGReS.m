function initGReS()
% List of folders to check
foldList = ["Code","ThirdPartyLibs","Utilities"];

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
  fprintf('Restoring default paths \n');
  restoredefaultpath;
  fprintf('Done Restoring default paths \n');
end

% set the root to the gres directory
gres_root = fileparts(mfilename('fullpath'));
setappdata(0,'gres_root', gres_root);

% Add GReS paths
for k = 1:numel(foldList)
  addpath(genpath(foldList(k)))
end

end


% addpath(genpath("Tests_Develop"));
