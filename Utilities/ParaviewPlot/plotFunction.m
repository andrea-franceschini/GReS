function  plotFunction(mesh, foldName, funct, varargin)
% PLOTFUNCTION Summary of this function goes here
%   Detailed explanation goes here
outVTK = VTKOutput(mesh, foldName); % create VTK object
pointData.name = 'solution';
pointData.data = funct;
if isempty(varargin) || strcmpi(varargin{1},'node')
   if ~isempty(mesh.cells)
      outVTK.writeVTKFile(0, pointData, [], [], []);
   else
      outVTK.writeVTKFile(0, [], [], pointData, []);
   end
elseif strcmpi(varargin{1},'elem')
   if ~isempty(mesh.cells)
      outVTK.writeVTKFile(0, [], pointData, [], []);
   else
      outVTK.writeVTKFile(0, [], [], [], pointData);
   end
   outVTK.finalize()
end

