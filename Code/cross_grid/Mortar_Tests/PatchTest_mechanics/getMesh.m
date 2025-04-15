function meshOut = getMesh(dir,fIn,meshName,nx,ny,varargin)
% Read and modify a geo file by setting number of nodes in each dimension
dirName = dir;
geoFile = strcat(dir,'/',fIn); % name of the geo file
% ntip is the number of non mortar nodes before the lateral boundary


% Read the file
fileLines = {};
fid = fopen(geoFile, 'r');
if fid == -1
    error('Cannot open the file.');
end

tline = fgetl(fid);
while ischar(tline)
    fileLines{end+1, 1} = tline; % Store each line in a cell array
    tline = fgetl(fid);
end
fclose(fid);

% Step 2: Modify mesh size and .msh file name
for i = 1:numel(fileLines)
    strVec = regexp(fileLines{i,1}, '(\S+|\s+)', 'match');
    if ~isempty(strVec)
        switch strVec{1}
            case 'nX'
                strVec{end} = strcat(num2str(nx),';');
            case 'nY'
                strVec{end} =  strcat(num2str(ny),';');
            case 'nXSlave'
                assert(~isempty(varargin),'No input for tip elements');
                strVec{end} =  strcat(num2str(varargin{1}),';');
            case 'nTip'
                assert(~isempty(varargin),'No input for tip elements');
                strVec{end} =  strcat(num2str(varargin{2}),';');
        end
    end
    fileLines{i,1} = strjoin(strVec, '');
end

% Step 3: Write back to the file
fid = fopen(geoFile, 'w');
if fid == -1
    error('Cannot open the file for writing.');
end

for i = 1:length(fileLines)
    fprintf(fid, '%s\n', fileLines{i});
end
fclose(fid);


% run the .geo file to produce the .msh file
meshFile =strcat(dirName,'/',meshName,'.msh');
command = strcat("LD_LIBRARY_PATH=; gmsh -2 ",geoFile," -o ",meshFile);
status = system(command);

% Check if the command was successful
% if status ~= 0
%     error('Gmsh not ran successfully');
% else

% create the mesh object 
meshOut = Mesh();
% gmsh -2 example.geo -nopopup -o example.msh
meshOut.importGMSHmesh(meshFile);

end

