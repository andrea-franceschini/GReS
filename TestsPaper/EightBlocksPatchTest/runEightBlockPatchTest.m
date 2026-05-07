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
grid = cell(8,1);

for i = 1:numel(meshList)
  meshList{i} = Mesh();
  meshList{i}.importMesh(fName + num2str(i) + ".vtk");
  msh = meshList{i};
  elems = Elements(msh,3);
  grid{i} = struct('topology',msh,'cells',elems,'faces',[]);
  % msh = meshList{i};
  % plotFunction(msh,"mesh_"+num2str(i),ones(msh.nNodes,1))
end

simparams = SimulationParameters('simParam.xml');
mat = Materials('materials.xml');

% write BC files
bcs = setBCfiles(grid);

domainFile = 'domains.xml';
interfFile = 'interfaces_nodal.xml';


for i = 1:8
  printUtils = OutState(meshList{i},"folderName",strcat("Output/mesh",num2str(i)),"timeList",1,...
    "writeVtk",1);
  domains(i,1) = Discretizer('grid',grid{i},...
    'materials',mat,...
    'boundaries',bcs{i},...
    'outstate',printUtils);
  domains(i,1).addPhysicsSolver('solver.xml');
end

interfaces = buildInterfaces(interfFile,domains);

solver = MultidomainFCSolver(simparams,domains,interfaces);


solver.NonLinearLoop();
solver.finalizeOutput();
