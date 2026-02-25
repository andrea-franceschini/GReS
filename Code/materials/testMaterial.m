clc 
clear

mat1 = Materials();
mat1.addSolid('name',"rock",'cellTags',1);

mat1.addConstitutiveLaw("rock","Elastic",'youngModulus',5e3,'poissonRatio',0.25);
mat1.addFluid('dynamicViscosity',1e-3);
mat1.addPorousRock("rock","specificWeight",21.0,"permeability",1e-12);

%%

mat = Materials("Dev/materials.xml");

