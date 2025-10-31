outFileTop = "topBlock";
outFileBottom = "bottomBlock";


domainFile = 'Domains/domains.xml';
interfFile = 'Domains/interface.xml';
domains = buildModel(domainFile); 
% set verbosity
domains(1).simparams.setVerbosity(0);
domains(2).simparams.setVerbosity(0);

[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);

solver = ActiveSetContactSolver(domains,interfaces,5);

%solver.simParameters.setBackstepSkipFlag(1);

solver.NonLinearLoop();
solver.finalizeOutput();

% get tangential gap
gt = interfaces{1}.tangentialGap.curr;
anGt = 0.1*sqrt(2);
tol = 1e-6;
assert(all(abs(gt - anGt)<tol),"Analytical solution is not matched")




% 

