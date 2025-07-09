function modelStruct = buildModelStruct_new(fileName,simParam)
% % build a structure array containing istances of all classes for different
% noncofnorming domains in the simulation
outStruct = readstruct(fileName,FileType="xml");
modelStruct = [];
for i = 1:numel(outStruct.Domain)
  modelStruct = [modelStruct; getModelStruct(outStruct.Domain(i),simParam)];
end
end


function modStr = getModelStruct(str,simParam)

  name = str.Name;
  model = ModelType(str.ModelType);
  topology = Mesh();
  topology.importMesh(char(str.Geometry));
  if isfield(str,"Material")
    material = Materials(model,char(str.Material));
  else
    material = [];
  end

  % dof manager
  if ~isfield(str,'DoFManager')
    dof = DoFManager(topology,model);
  else
    dof = DoFManager(topology,model,str.DoFManager);
  end
  
  gNPoints = str.Gauss;
  
  if ~ismissing(gNPoints)
    elems = Elements(topology,gNPoints);
  else
    elems = Elements(topology);
  end
  faces = Faces(model,topology);
  grid = struct('topology',topology,'cells',elems,'faces',faces);
 
  % boundary conditions
  nBC = str.BoundaryConditions.countAttribute;
  bcList = strings(1,nBC);
  for i = 1:nBC
    bcList(i) = str.BoundaryConditions.File(i);
  end
  bc = Boundaries(bcList,model,grid);

  % output manager
  printUtils = OutState(model,topology,str.OutState,"folderName",name,"writeVtk",false);

  linSyst = Discretizer(model,simParam,dof,grid,material);

  linSyst.initState();

  modStr = struct( ...
    'DomainName',         name, ...
    'ModelType',          model, ...
    'Grid',               grid, ...
    'Material',           material, ...
    'DoFManager',         dof, ...
    'BoundaryConditions', bc, ...
    'OutState',           printUtils, ...
    'Discretizer',        linSyst ...
);
end

