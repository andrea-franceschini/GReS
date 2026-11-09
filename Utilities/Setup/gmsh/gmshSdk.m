function [sdk, ready] = gmshSdk()
% Describe the pinned SDK and check the files needed to link and load it.
version = "4.15.2";
abiVersion = regexprep(version, '\.[^.]+$', '');
switch computer
  case 'PCWIN64'
    folder = "Windows";
    platform = "Windows64";
    extension = ".zip";
    linkName = "gmsh.dll.lib";
    runtimeName = "gmsh-" + abiVersion + ".dll";
  case 'GLNXA64'
    folder = "Linux";
    platform = "nox-Linux64";
    extension = ".tgz";
    linkName = "libgmsh.so";
    runtimeName = "libgmsh.so." + abiVersion;
  case 'MACA64'
    folder = "MacOSX";
    platform = "MacOSARM";
    extension = ".tgz";
    linkName = "libgmsh.dylib";
    runtimeName = "libgmsh." + abiVersion + ".dylib";
  otherwise
    sdk = [];
    ready = false;
    return
end

% The archive's top-level directory records the version and architecture.
sdk.name = "gmsh-" + version + "-" + platform + "-sdk";
sdk.archive = sdk.name + extension;
sdk.url = "https://gmsh.info/bin/" + folder + "/" + sdk.archive;
sdk.parent = fullfile(gres_root, 'ThirdPartyLibs', 'Gmsh');
sdk.directory = fullfile(sdk.parent, sdk.name);
sdk.includeDir = fullfile(sdk.directory, 'include');
sdk.libDir = fullfile(sdk.directory, 'lib');
sdk.linkFile = fullfile(sdk.libDir, linkName);
sdk.requiredFiles = ["include/gmshc.h", "lib/" + linkName, "lib/" + runtimeName];
ready = all(isfile(fullfile(sdk.directory, sdk.requiredFiles)));
end
