mat = Materials();
mat.addSolid('name',"rock",'cellTags',1);

mat.addConstitutiveLaw(1,"Elastic",'youngModulus',5e3,'poissonRatio',0.25);

%%

mat = Materials("materials.xml");
