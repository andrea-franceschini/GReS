function  plotFunction(mesh, foldName, funct, name, varargin)
% shortcut function to produce VTK output
% array solutions can be handled
% deal inpts
t = 0;
varType = [];
outVTK = [];
for j = 1:numel(varargin)
   switch varargin{j}
      case 'time'
         t = varargin{j+1};
      case 'type'
         varType = varargin{j+1};
      case 'VTK'
         outVTK = varargin{j+1};
   end
end
if isempty(outVTK)
   outVTK = VTKOutput(mesh,foldName);
end
ncomp = size(funct,2);
assert(numel(name)==ncomp,'Wrong array size for input name');

pointData = repmat(struct('name',[],'data',[]),ncomp,1);
for i = 1:size(funct,2)
    pointData(i).name = convertStringsToChars(name(i));
    pointData(i).data = funct(:,i);
end
if isempty(varargin) || strcmpi(varType,'node')
   if ~isempty(mesh.cells)
      outVTK.writeVTKFile(t, pointData, [], [], []);
   else
      outVTK.writeVTKFile(t, [], [], pointData, []);
   end
elseif strcmpi(varType,'elem')
   if ~isempty(mesh.cells)
      outVTK.writeVTKFile(t, [], pointData, [], []);
   else
      outVTK.writeVTKFile(t, [], [], [], pointData);
   end
end
   outVTK.finalize()
end