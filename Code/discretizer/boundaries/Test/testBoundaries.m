fileName = 'boundaryConditions.xml';

model = ModelType(["SinglePhaseFlow_FVTPFA","Poromechanics_FEM"]);

mesh = Mesh();
mesh.importMesh("Column_hexa.msh");
elems = Elements(mesh,2);
faces = Faces(model,mesh);

grid = struct("topology",mesh,"cells",elems,"faces",faces);

bc = BoundariesNew(fileName,model,grid);

%% test some bc methods
bcList = bc.db.keys; 

v1 = getVals(bc,"TopLoad",0.5);
v2 = getVals(bc,"TopLoad",1.5);
v3 = getVals(bc,"TopLoad",8);


v4 =  getVals(bc,"LatFixedX",100);

e1 = getEntities(bc,"NoFlowTop");
e2 = getEntities(bc,"TopLoad");
e3 = getEntities(bc,"LatFixedX");


% for i = 1:numel(bcList)
%   id = bcList{i};
%   cond = bc.getCond(id);
%   ents = bc.getEntities(id);
%   entInf = bc.getEntitiesInfluence(id);
% end