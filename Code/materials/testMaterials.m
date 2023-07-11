close all;
clear;

% Setting the input file name
fileName = 'materialsList.dat';

% Creation and intialization of an object of "Materials"
mat = Materials(fileName);

% Calling the "getMaterial" function from the class "Materials" 
D1 = mat.getMaterial(1).ConstLaw.getStiffnessMatrix()
K1 = mat.getMaterial(1).PorousRock.getPermMatrix()
por1 = mat.getMaterial(1).PorousRock.getPorosity()
swr1 = mat.getMaterial(1).PorousRock.getWaterResSat()
%
% D2 = mat.getMaterial(2).ConstLaw.getStiffnessMatrix();
K2 = mat.getMaterial(2).PorousRock.getPermMatrix()
por2 = mat.getMaterial(2).PorousRock.getPorosity()
swr2 = mat.getMaterial(2).PorousRock.getWaterResSat()
%
gWater = mat.getMaterial(3).getFluidSpecWeight()
betaWater = mat.getMaterial(3).getFluidCompressibility()