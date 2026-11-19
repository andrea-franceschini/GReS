function available = isAvailable()
%ISAVAILABLE Check for the Gmsh MEX built for this MATLAB platform.
gateway = fullfile(fileparts(mfilename('fullpath')), 'private', ...
  ['gresGmshGateway.' mexext]);
available = exist(gateway, 'file') == 3;
end
