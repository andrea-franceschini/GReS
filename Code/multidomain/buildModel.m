function domains = buildModel(fileName)
% % build a structure array containing istances of all classes for different
% noncofnorming domains in the simulation
outStruct = readstruct(fileName,AttributeSuffix="");
if isfield(outStruct,"Domains")
  outStruct = outStruct.Domains;
end

if isfield(outStruct,"Domain")
 outStruct = outStruct.Domain;
end

nD = numel(outStruct);
domains = cell(nD,1);

for i = 1:nD
  domains{i} = defineDomain(outStruct(i));
end

domains = cell2mat(domains);

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


  gNPoints = getXMLData(struc,0,'Gauss');

  if gNPoints~=0
    elems = Elements(topology,gNPoints);
  else
    elems = Elements(topology);
  end

  faces = Faces(model,topology);
  grid = struct('topology',topology,'cells',elems,'faces',faces);

  % boundary conditions
  if isfield(struc,'BoundaryConditions')
    bcFile = getXMLData(struc,[],"BoundaryConditions");
    bc = Boundaries(bcFile,model,grid);
  else
    bc = [];
  end

  % output manager
  if isfield(struc,"OutState")
    printFlag = getXMLData(struc.OutState,0,'print');
    matFileFlag = getXMLData(struc.OutState,0,'matFile');
    outFile = getXMLData(struc.OutState,[],'outFile');
    printUtils = OutState(model,topology,outFile,...
      "folderName",name,"writeVtk",printFlag,"flagMatFile",matFileFlag);
  end

  modStr = Discretizer('ModelType',model,...
                      'SimulationParameters',simparam,...
                      'DoFManager',dof,...
                      'Boundaries',bc,...
                      'OutState',printUtils,...
                      'Materials',material,...
                      'Grid',grid);


end

