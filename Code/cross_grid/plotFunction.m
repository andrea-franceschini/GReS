function  plotFunction(mesh, foldName, funct)
%PLOTFUNCTION Summary of this function goes here
%   Detailed explanation goes here
outVTK = VTKOutput(mesh, foldName); % create VTK object
pointData.name = 'solution';
pointData.data = funct;
if any(mesh.coordinates(:,3)~=0)
    outVTK.writeVTKFile(0, pointData, [], [], []);
else
outVTK.writeVTKFile(0, [], [], pointData, []);
end
outVTK.finalize()
end

