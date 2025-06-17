function initGReS()


% Remove existing paths
restoredefaultpath;

% set the root to the gres directory
gres_root = fileparts(mfilename('fullpath'));
setappdata(0,'gres_root', gres_root);

% Init GRES
addpath(genpath("Code"));
addpath(genpath("ThirdPartyLibs"));
addpath(genpath("Utilities"));
%addpath(genpath("Tests_Develop"));
end
