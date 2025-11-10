clear
close all
clc


 
scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);
%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('Eight blocks patch test \n')
fprintf('___________________\n\n')


fName = "Mesh/cube";

meshList = cell(8,1);
for i = 1:numel(meshList)
  meshList{i} = Mesh();
  meshList{i}.importMesh(fName + num2str(i) + ".vtk");
  msh = meshList{i};
  plotFunction(msh,"mesh_"+num2str(i),ones(msh.nNodes,1))
end


% write BC files
setBCfiles(meshList);

domainFile = 'domains.xml';
interfFile = 'interfaces.xml';
domains = buildModel(domainFile); 

[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);

solver = MultidomainFCSolver(domains,interfaces);


solver.NonLinearLoop();
solver.finalizeOutput();
