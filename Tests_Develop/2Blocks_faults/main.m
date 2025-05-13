clear
close all

% fault material parameters
coes = 0;
phi = 25; % degrees

multType = 'P0';

fprintf('\n ________________ \n')
fprintf('%s multipliers \n',multType)
% Get the full path of the script
scriptFullPath = mfilename('fullpath');
% Extract the directory of the script
scriptDir = fileparts(scriptFullPath);
% Change to the script's directory
cd(scriptDir);
%str_mod_ver = '2km'; % '2km' or '3km'

fprintf('2blocks model \n')
fprintf('___________________\n\n')
[leftMesh,rightMesh] = deal(Mesh(),Mesh());

simulTag = 'hexa';
meshLeftFile = 'Mesh/LeftBlock_hexa.msh';
meshRightFile = 'Mesh/RightBlock_hexa.msh';
domainFile = 'Domains/domains_hexa_P0.xml';

leftMesh.importGMSHmesh(meshLeftFile);
rightMesh.importGMSHmesh(meshRightFile);
% 
% plotFunction(leftMesh,'out_L',zeros(leftMesh.nNodes,1),"sol");
% plotFunction(rightMesh,'out_R',zeros(rightMesh.nNodes,1),"sol");

% write BC files
setBCfiles(leftMesh,rightMesh);

interfFile = 'interface.xml';
simParam = SimulationParameters('simParam.dat');
modelStruct = buildModelStruct_new(domainFile,simParam);
interfaces = MeshGlue.buildInterfaceStruct(interfFile,modelStruct);


% validate mortar operator construction
% f = @(y,z) sin(0.5*y)+cos(0.5*z);
% cs = mG.interfaces.mortar.intSlave.coordinates;
% cm = mG.interfaces.mortar.intMaster.coordinates;
% fSAnal = f(cs(:,2),cs(:,3));
% fMaster = f(cm(:,2),cm(:,3));
% fSInterp = mG.interfaces.InterpOperator*fMaster;
% e = fSInterp-fSAnal;
% e = sqrt(sum(e.^2));

solverDual = NonLinearSolverFaults(simParam,mG,coes,phi,multType);

solverDual.models(1).OutState.finalize();
solverDual.models(2).OutState.finalize();
% measure accuracy
% errMult = sqrt(sum((abs(solver.currMultipliers(1:3:end))-2).^2));
%% POST PROCESSING

% % Define the folder path
% folderPath = 'OUT';
% if strcmp(multType,'dual')
%    % List all files in the folder with the specified pattern
%    filePattern = fullfile(folderPath, '*dual*');
% elseif strcmp(multType,'standard')
%    filePattern = fullfile(folderPath, '*standard*');
% elseif strcmp(multType,'P0')
%    filePattern = fullfile(folderPath, '*P0*');
% end
% 
% % Get all entries matching the pattern
% entriesToDelete = dir(filePattern);
% 
% % Loop through and delete each entry
% for k = 1:length(entriesToDelete)
%     entryName = entriesToDelete(k).name;
%     entryPath = fullfile(folderPath, entryName);
%     % Check if it's a file or directory
%     if isfile(entryPath)
%         delete(entryPath); % Delete the file
%     elseif isfolder(entryPath)
%         rmdir(entryPath, 's'); % Delete the directory and its contents
%     end
% end
% 
% 
% movefile Left* OUT
% movefile Right* OUT
% movefile Fault* OUT 


%%

% multType = 'standard';
% 
% fprintf('\n ________________ \n')
% fprintf('%s multipliers \n',multType)
% % Get the full path of the script
% scriptFullPath = mfilename('fullpath');
% % Extract the directory of the script
% scriptDir = fileparts(scriptFullPath);
% % Change to the script's directory
% cd(scriptDir);
% %str_mod_ver = '2km'; % '2km' or '3km'
% 
% fprintf('2blocks model \n')
% fprintf('___________________\n\n')
% [leftMesh,rightMesh] = deal(Mesh(),Mesh());
% 
% simulTag = 'hexa';
% switch simulTag
%    case 'tetra'
%       meshLeftFile = 'Mesh/LeftBlock_tetra.msh';
%       meshRightFile = 'Mesh/RightBlock_tetra.msh';
%       switch multType
%          case 'standard'
%             domainFile = 'Domains/domains_tetra_standard.dat';
%          case 'dual'
%             domainFile = 'Domains/domains_tetra_dual.dat';
%       end
%    case 'hexa'
%       meshLeftFile = 'Mesh/LeftBlock_hexa.msh';
%       meshRightFile = 'Mesh/RightBlock_hexa.msh';
%       switch multType
%          case 'standard'
%             domainFile = 'Domains/domains_hexa_standard.dat';
%          case 'dual'
%             domainFile = 'Domains/domains_hexa_dual.dat';
%       end
% end
% 
% leftMesh.importGMSHmesh(meshLeftFile);
% rightMesh.importGMSHmesh(meshRightFile);
% % 
% % plotFunction(leftMesh,'out_L',zeros(leftMesh.nNodes,1),"sol");
% % plotFunction(rightMesh,'out_R',zeros(rightMesh.nNodes,1),"sol");
% 
% % write BC files
% setBCfiles(leftMesh,rightMesh);
% 
% switch multType
%    case 'dual'
%       interfFile = 'interface_dual.dat';
%    case 'standard'
%       interfFile = 'interface_standard.dat';
% end
% simParam = SimulationParameters('simParam.dat');
% model = buildModelStruct(domainFile,simParam);
% %scale domain size to millimiters
% mG = MeshGlue(model,interfFile);
% 
% 
% % validate mortar operator construction
% % f = @(y,z) sin(0.5*y)+cos(0.5*z);
% % cs = mG.interfaces.mortar.intSlave.coordinates;
% % cm = mG.interfaces.mortar.intMaster.coordinates;
% % fSAnal = f(cs(:,2),cs(:,3));
% % fMaster = f(cm(:,2),cm(:,3));
% % fSInterp = mG.interfaces.InterpOperator*fMaster;
% % e = fSInterp-fSAnal;
% % e = sqrt(sum(e.^2));
% 
% solverStandard = NonLinearSolverFaults(simParam,mG,coes,phi);
% 
% solverStandard.models(1).OutState.finalize();
% solverStandard.models(2).OutState.finalize();
% measure accuracy
% errMult = sqrt(sum((abs(solver.currMultipliers(1:3:end))-2).^2));
%% POST PROCESSING

%Define the folder path
folderPath = 'OUT';
if strcmp(multType,'dual')
   %List all files in the folder with the specified pattern
   filePattern = fullfile(folderPath, '*dual*');
elseif strcmp(multType,'standard')
   filePattern = fullfile(folderPath, '*standard*');
   elseif strcmp(multType,'P0')
   filePattern = fullfile(folderPath, '*P0*');
end

% Get all entries matching the pattern
entriesToDelete = dir(filePattern);

% Loop through and delete each entry
for k = 1:length(entriesToDelete)
    entryName = entriesToDelete(k).name;
    entryPath = fullfile(folderPath, entryName);
    % Check if it's a file or directory
    if isfile(entryPath)
        delete(entryPath); % Delete the file
    elseif isfolder(entryPath)
        rmdir(entryPath, 's'); % Delete the directory and its contents
    end
end

movefile Left* OUT
movefile Right* OUT
movefile Fault* OUT 

%% plot profiles of multipliers along vertical axis (avoid opening paraview)
plotStep(solverDual.results,6);




% 

