model = ModelManager();
model.createModel('terzaghi.xml');

% option 2
% model.addDomain('grid',grid,...
%                 'materials',mat,...
%                 'boundaryConditions',bcs, ...
%                 'printUtils',printUtils);
% 
% model.addDomain('grid',grid,...
%                 'materials',mat,...
%                 'boundaryConditions',bcs, ...
%                 'printUtils',printUtils);
% 
% model.getDomain(1).addPhysicsSolver("Poromechanics",'targetRegions',1);
% model.getDomain(2).addPhysicsSolver("Biot",'targetRegions',[1 2]);
% 
% model.solutionStrategy =  

model.runProblem();