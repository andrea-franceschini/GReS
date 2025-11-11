function domains = buildModel(fileName)
% % build a structure array containing istances of all classes for different
% noncofnorming domains in the simulation
outStruct = readstruct(fileName,FileType="xml");
domains = [];
for i = 1:numel(outStruct.Domain)
  domains = [domains; defineDomain(outStruct.Domain(i))];
end
end


function modStr = defineDomain(struc)

  name = struc.Name;
  model = ModelType(strsplit(struc.ModelType));
  topology = Mesh();
  topology.importMesh(char(struc.Geometry));
  if isfield(struc,"Material")
    material = Materials(char(struc.Material));
  else
    material = [];
  end

  % dof manager
  if ~isfield(struc,'DoFManager')
    dof = DoFManager(topology,model);
  else
    dof = DoFManager(topology,model,struc.DoFManager);
  end

  if isfield(struc,"SimParameters")
    simparam = SimulationParameters(struc.SimParameters);
  end


  gNPoints = getFieldAttribute(struc,'Gauss',[]);

  if ~isempty(gNPoints)
    elems = Elements(topology,gNPoints);
  else
    elems = Elements(topology);
  end

  faces = Faces(model,topology);
  grid = struct('topology',topology,'cells',elems,'faces',faces);

  % boundary conditions
  if isfield(struc,'BoundaryConditions')
    nBC = struc.BoundaryConditions.countAttribute;
    bcList = strings(1,nBC);
    for i = 1:nBC
      bcList(i) = struc.BoundaryConditions.File(i);
    end
    bc = Boundaries(bcList,model,grid);
  else
    bc = [];
  end

  % output manager
  if isfield(struc,"OutState")
    printFlag = getFieldAttribute(struc.OutState,'printAttribute',0);
    matFileFlag = getFieldAttribute(struc.OutState,'matFileAttribute',0);
    printUtils = OutState(model,topology,struc.OutState.outFileAttribute,...
      "folderName",name,"writeVtk",printFlag,"flagMatFile",matFileFlag);
  else
    printUtils = [];
  end

  modStr = Discretizer('ModelType',model,...
                      'SimulationParameters',simparam,...
                      'DoFManager',dof,...
                      'Boundaries',bc,...
                      'OutState',printUtils,...
                      'Materials',material,...
                      'Grid',grid);


end

function out = getFieldAttribute(str,attr,default)
if isfield(str,attr)
  out = str.(attr);
else
  if nargin > 2
    out = default;
  else
    error('Missingrequired attribute %s in domain file \n', attr)
  end
end


end