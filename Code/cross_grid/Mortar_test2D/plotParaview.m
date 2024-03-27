function  plotParaview(mesh, foldName, u, dir)
%PLOTFUNCTION Summary of this function goes here
%   Detailed explanation goes here
outVTK = VTKOutput(mesh, foldName); % create VTK object

if strcmp(dir, 'all')
pointData(1,1).name = 'ux';
pointData(1,1).data = u(1:2:end)';
pointData(2,1).name = 'uy';
pointData(2,1).data = u(2:2:end)';
elseif strcmp(dir, 'x')
    pointData(1,1).name = 'ux';
    pointData(1,1).data = u';
elseif strcmp(dir, 'y')
    pointData(1,1).name = 'uy';
    pointData(1,1).data = u';    
end
if any(mesh.coordinates(:,3)~=0)
    outVTK.writeVTKFile(0, pointData, [], [], []);
else
outVTK.writeVTKFile(0, [], [], pointData, []);
end
outVTK.finalize()
end

