function  plotMesh(mesh, foldName, funct, varargin)
% PLOTFUNCTION Summary of this function goes here
%   Detailed explanation goes here
outVTK = VTKOutput(mesh, foldName); % create VTK object
outVTK.writeVTKFile(0, [], [], [], []);
outVTK.finalize()
end

