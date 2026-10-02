clear;
clc;

% Empties all Output/Outputs folders, preserving the folders themselves.
% Deletes all .mat files larger than the specified threshold.

sizeLimitMB = 0.25;                 % MB = 1,000,000 bytes
sizeLimitBytes = sizeLimitMB * 1e6;
repositoryRoot = fileparts(fileparts(mfilename('fullpath')));

entries = dir(fullfile(repositoryRoot, '**', '*'));

% Find Output/Outputs folders, case insensitive.
isOutput = [entries.isdir] & ...
    (strcmpi({entries.name}, 'Output') | ...
     strcmpi({entries.name}, 'Outputs'));

outputEntries = entries(isOutput);
outputFolders = strings(numel(outputEntries), 1);
for k = 1:numel(outputEntries)
    outputFolders(k) = fullfile(outputEntries(k).folder, outputEntries(k).name);
end

% Retain only outermost Output folders to avoid duplicate deletions.
[~, order] = sort(strlength(outputFolders));
outputFolders = outputFolders(order);
selectedFolders = strings(0, 1);

for k = 1:numel(outputFolders)
    if ~any(startsWith(outputFolders(k), selectedFolders + filesep))
        selectedFolders(end + 1, 1) = outputFolders(k); %#ok<SAGROW>
    end
end

% Find large MAT files outside the folders already scheduled for cleaning.
isMat = ~[entries.isdir] & ...
    endsWith({entries.name}, '.mat', 'IgnoreCase', true) & ...
    [entries.bytes] > sizeLimitBytes;

matEntries = entries(isMat);
matFiles = strings(0, 1);

for k = 1:numel(matEntries)
    filePath = string(fullfile(matEntries(k).folder, matEntries(k).name));
    if ~any(startsWith(filePath, selectedFolders + filesep))
        matFiles(end + 1, 1) = filePath; %#ok<SAGROW>
    end
end

fprintf('Repository: %s\n', repositoryRoot);
fprintf('\nOutput folders whose contents will be deleted:\n');
disp(selectedFolders);

fprintf('Additional MAT files larger than %.3f MB:\n', sizeLimitMB);
disp(matFiles);

if isempty(selectedFolders) && isempty(matFiles)
    fprintf('Nothing to delete.\n');
    return
end

answer = '';
while ~ismember(answer, {'y', 'n'})
    answer = lower(strtrim(input( ...
        'Permanently delete the listed contents and files? [y/n]: ', 's')));
end

if strcmp(answer, 'n')
    fprintf('Cancelled. No files deleted.\n');
    return
end

% Empty Output folders while preserving each Output folder.
for k = 1:numel(selectedFolders)
    children = dir(selectedFolders(k));
    children = children(~ismember({children.name}, {'.', '..'}));

    for j = 1:numel(children)
        target = fullfile(children(j).folder, children(j).name);
        if children(j).isdir
            rmdir(target, 's');
        else
            delete(target);
        end
    end
end

% Delete the additional large MAT files.
for k = 1:numel(matFiles)
    if isfile(matFiles(k))
        delete(matFiles(k));
    end
end

fprintf('Repository cleanup completed.\n');
