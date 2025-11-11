fileName = gres_root + "/docs/inputExamples/boundaryConditions.xml";

model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);

mesh = Mesh();
mesh.importMesh("Column_hexa.msh");
elems = Elements(mesh,2);


grid = struct("topology",mesh,"cells",elems);

bc = BoundariesNew(fileName,model,grid);

%% test some bc methods
bcList = bc.db.keys; 

for i = 1:numel(bcList)
  id = bcList{i};
  cond = bc.getCond(id);
  ents = bc.getEntities(id);
  entInf = bc.getEntitiesInfluence(id);
end