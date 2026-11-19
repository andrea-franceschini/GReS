function compileGmsh(required)
%COMPILEGMSH Build the native gateway with the installed Gmsh SDK.
[sdk, ready] = gmshSdk();
if ~ready
  if required
    error('GReS:GmshSDKUnavailable', ...
      'Gmsh SDK unavailable for %s. Run setupGmsh to enable Gmsh support.', computer);
  end
  fprintf('Gmsh SDK unavailable for %s. Skipping Gmsh interface.\n', computer);
  fprintf('Run setupGmsh to enable Gmsh support.\n');
  return
end

mexDir = fileparts(mfilename('fullpath'));
gatewayDir = fullfile(mexDir, '+gresgmsh', 'private');
if ~isfolder(gatewayDir)
  mkdir(gatewayDir);
end
args = ["-R2017b", "-O", "-outdir", string(gatewayDir), "-output", "gresGmshGateway", ...
  "-I" + sdk.includeDir, fullfile(mexDir,'gmshMex.cpp'), sdk.linkFile];
if ~ispc
  args(end+1) = "CXXFLAGS=$CXXFLAGS -std=c++14";
  args(end+1) = 'LDFLAGS=$LDFLAGS -Wl,-rpath,"' + sdk.libDir + '"';
end
fprintf('Compiling Gmsh MEX with %s\n', sdk.name);
args = cellstr(args);
mex(args{:});
addGmshRuntimePath(sdk);
end
