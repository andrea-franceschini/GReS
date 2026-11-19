function addGmshRuntimePath(sdk)
% Windows locates the SDK DLL through MATLAB's process-local PATH.
if ispc
  currentPath = getenv('PATH');
  if ~any(strcmpi(split(string(currentPath), pathsep), sdk.libDir))
    setenv('PATH', sdk.libDir + pathsep + currentPath);
  end
end
end
