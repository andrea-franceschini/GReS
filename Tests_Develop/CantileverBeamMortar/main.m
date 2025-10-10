clear
close all
clc


scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);
%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('Cantilever beam validation for GEOS \n')
fprintf('___________________\n\n')

mesh = Mesh();

mesh.importMesh("Mesh/beam.vtk");

writeBCfiles('BCs/fix','SurfBC','Dir',{'Poromechanics','x','y','z'},'1fix',0,0,mesh,3); % left block lateral fix
writeBCfiles('BCs/topLoad','SurfBC','Neu',{'Poromechanics','z'},'0load',0,-1e-2,mesh,4); % left block bottom fix


domainFile = 'Domains/domains.xml';
interfFile = 'Domains/interface.xml';
domains = buildModel(domainFile); 
% set verbosity 
[interfaces,domains] = Mortar.buildInterfaces(interfFile,domains);


interfaces{1}.solvers(2).simparams.setVerbosity(2);

solver = MultidomainFCSolver(domains,interfaces);


%solver.simParameters.setBackstepSkipFlag(1);

solver.NonLinearLoop();
solver.finalizeOutput();

%% plot profiles of multipliers along vertical axis (avoid opening paraview)
%plotStep(solver.results,2);


% matrixGEOS = load("systemMatrix.mtx");
% 
% nMult = 100;
% 
% r = matrixGEOS(2:end,1);
% c = matrixGEOS(2:end,2);
% v = matrixGEOS(2:end,3);
% 
% mat = sparse(r,c,v);
% figure(1)
% %title("traction_bubble")
% spy(mat)
% 
% nDof = size(mat,1);
% 
% tDof = nDof-nMult*3*3;
% 
% bDof = nDof-nMult*2*3;
% 
% % get the Atb block (probably the wrong one)
% Abt_geos = mat(tDof+1:tDof+nMult*3, bDof+1:bDof+nMult*3);
% Abb_geos = mat(bDof+1:bDof+nMult*3,bDof+1:bDof+nMult*3);
% 






% 

