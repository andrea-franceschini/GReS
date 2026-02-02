function  plotMesh(mesh, foldName)
% PLOTFUNCTION Summary of this function goes here
%   Detailed explanation goes here
outVTK = VTKOutput(mesh, foldName); % create VTK object
pointData.name = 'field';
pointData.data = zeros(mesh.nNodes,1);
outVTK.writeVTKFile(0, pointData, [], [], []);
outVTK.finalize()
end

