% testing the DruckerPrager driver for validation

% input parameter follow GEOS validation study



law = DruckerPrager('param.xml');


driver = TriaxialDriver('nStep',200,'constLaw',law,'outFile',"DPtest.txt");

t = 0:5;
axStrain = -[0;4;2;5;3;6]*1e-3;
radStress = -10e6*ones(6,1);
driver.setFunction('axialStrain',t,axStrain);
driver.setFunction('radialStress',t,radStress);

driver.setParameters('initialStress',-10e6);

out = driver.launch;

validateDruckerPragerAnalytical(law,out,'StressUnit','Pa');

