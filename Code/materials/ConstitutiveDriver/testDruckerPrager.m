% testing the DruckerPrager driver for validation

% input parameter follow GEOS validation study



law = DruckerPrager('param.xml');


driver = TriaxialDriver('nStep',200,'constLaw',law,'outFile',"DPtest.txt");

driver.setFunction('');
driver.setFunction('');

