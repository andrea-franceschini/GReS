function  plotFunction(mesh, foldName, time, funct, varargin)
% PLOTFUNCTION Summary of this function goes here
%   Detailed explanation goes here

out = OutState('folderName',foldName,'timeList',time);
out.prepareOutputFolders();

[~, vtuName, ~] = fileparts(foldName);

toc = out.vtkFile.getDocumentElement;
toc.setAttribute('type', 'vtkMultiBlockDataSet');
toc.setAttribute('version', '1.0');
blocks = out.vtkFile.createElement('vtkMultiBlockDataSet');


% inner block with vtu dataset to append to vtm file
block = out.vtkFile.createElement('Block');

var.name = 'solution';
var.data = funct;
% outVTK.writeVTKFile(0, [], [], [], []);
if isempty(varargin) || strcmpi(varargin{1},'node')
  if ~isempty(mesh.cells)
    writeVTKfile(out,block,vtuName,mesh,time,var,[],[],[]);
  else
    writeVTKfile(out,block,vtuName,mesh,time,[],[],var,[]);
  end
elseif strcmpi(varargin{1},'elem')
  if ~isempty(mesh.cells)
    writeVTKfile(out,block,vtuName,mesh,time,[],var,[],[]);
  else
    writeVTKfile(out,block,vtuName,mesh,time,[],[],[],var);
  end
end

blocks.appendChild(block);
toc.appendChild(blocks);

out.writeVTMFile();

end

