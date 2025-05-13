function getPatchMesh(fIn,meshName,Nx,Ny)
% Read and modify a geo file by setting number of nodes in each dimension
geoFile = fIn; % name of the geo file
fInsplt = strsplit(fIn,'/');
dirName = fInsplt{1};

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
            case 'NX'
                strVec{end} = strcat(num2str(Nx),';');
           case 'NY'
                strVec{end} = strcat(num2str(Ny),';');
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
if isunix
  command = strcat("LD_LIBRARY_PATH=; gmsh -2 ",geoFile," -o ",meshFile);
elseif ispc
  command = strcat("C:\Users\dmoretto\Downloads\gmsh-4.13.1-Windows64\gmsh-4.13.1-Windows64\gmsh.exe -2 ",geoFile," -o ",meshFile);
end

status = system(command);

end
