clear
clc

% input parameter follow GEOS validation study

law = BiElastic('paramBielastic.xml');

driver = TriaxialDriver('nStep',300,'constLaw',law,'outFile',"BiElasticTest.txt");

axStress = -[0;10;3;20;6;30]*1e6;
radStress= axStress;
t = 0:(length(axStress)-1);

driver.setFunction('axialStress',t,axStress);
driver.setFunction('radialStress',t,radStress);
driver.setParameters('initialStress',0.0);

out = driver.launch;

fprintf('Done \n')

%validateDruckerPragerAnalytical(law,out,'StressUnit','Pa');


%% plotting
stressAx = out(:,driver.SIG0);
strainAx = out(:,driver.EPS0);

plot(-stressAx,-strainAx,'k-o')
xlabel('axial stress')
ylabel('axial strain')




