function writePlotCollection(grids, foldName, time, values, location, names)
% Use the simulation output layout for a single snapshot.
out = OutState('outputFile',foldName,'printTimes',time,'vtkFormat',"ascii");
out.prepareOutputFolders();

document = out.vtkFile;
root = document.getDocumentElement();
root.setAttribute('type','vtkMultiBlockDataSet');
root.setAttribute('version','1.0');
blocks = document.createElement('vtkMultiBlockDataSet');
root.appendChild(blocks);

for k = 1:numel(grids)
  grid = grids{k};
  block = document.createElement('Block');
  block.setAttribute('index',num2str(k-1));
  block.setAttribute('name',char(names(k)));
  blocks.appendChild(block);

  % Use the grid's populated topology. FractureGrid stores triangles in cells.
  data = {[],[],[],[]};
  slot = 1 + (location == "elem") + 2*(grid.cells.num == 0);
  data{slot} = struct('name','solution','data',values{k});
  writeVTKfile(out,block,sprintf('grid_%04d',k),grid,time,data{:});
  dataset = block.getLastChild();
  dataset.setAttribute('index','0');
end

out.writeVTMFile();
out.savePvd(time);
end
