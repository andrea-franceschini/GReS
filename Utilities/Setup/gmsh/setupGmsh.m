function setupGmsh(options)
%SETUPGMSH Install the optional, pinned Gmsh SDK for this MATLAB platform.
%   setupGmsh('force', true) downloads and replaces an existing installation.
arguments
  options.force (1,1) logical = false
end

[sdk, ready] = gmshSdk();
if isempty(sdk)
  error('GReS:GmshUnsupportedPlatform', 'No Gmsh SDK configured for %s.', computer);
end
if ready && ~options.force
  fprintf('%s is already installed.\n', sdk.name);
  return
end

if ~isfolder(sdk.parent)
  mkdir(sdk.parent);
end
staging = tempname(sdk.parent);
mkdir(staging);
cleanup = onCleanup(@() rmdir(staging, 's'));
archive = fullfile(staging, sdk.archive);
try
  websave(archive, sdk.url, weboptions('Timeout', 120));
catch cause
  error('GReS:GmshDownloadFailed', 'Could not download %s.\n%s', sdk.url, cause.message);
end

if ispc
  unzip(archive, staging);
else
  % Preserve the relative shared-library symlinks in the Unix SDKs.
  shellQuote = @(value) "'" + replace(value, "'", "'\''") + "'";
  [status, message] = system("tar -xzf " + shellQuote(archive) + ...
    " -C " + shellQuote(staging));
  if status ~= 0
    error('GReS:GmshExtractFailed', 'Could not extract %s.\n%s', sdk.archive, message);
  end
end
extracted = fullfile(staging, sdk.name);
if ~all(isfile(fullfile(extracted, sdk.requiredFiles)))
  error('GReS:GmshInvalidSDK', 'Expected Gmsh headers and libraries in %s.', sdk.url);
end

backup = tempname(sdk.parent);
hadInstallation = isfolder(sdk.directory);
if hadInstallation
  movefile(sdk.directory, backup);
end
[installed, message] = movefile(extracted, sdk.directory);
if ~installed
  if hadInstallation
    movefile(backup, sdk.directory);
  end
  error('GReS:GmshInstallFailed', 'Could not install %s.\n%s', sdk.name, message);
end
if hadInstallation
  rmdir(backup, 's');
end
fprintf('Installed %s.\n', sdk.name);
end
