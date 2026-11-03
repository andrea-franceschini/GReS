function stats = countLines()
% countLines - Count lines of code in GReS and ChronosLab without external tools.
%
% Computes comment and code lines (excluding blank lines) separately for:
%   - GReS (Main repository)
%   - ChronosLab (submodule)
% Displayed in a single unified table per repository with:
%   - MATLAB (.m)
%   - C / C++ (sources & headers)
%   - MEX Gateways (entry points)
%   - TOTAL (Sum of all categories)
%
% Percentages are displayed inside parentheses against total count (comment + code).
% MRST is excluded.

    gresRoot = gres_root();

    fprintf('\n');
    s1 = analyzeRepo('GReS (Main Repository)', gresRoot);
    s2 = analyzeRepo('ChronosLab', fullfile(gresRoot, 'ThirdPartyLibs', 'ChronosLab'));

    if nargout > 0
        stats = struct('GReS', s1, 'ChronosLab', s2);
    end
end

function repoStats = analyzeRepo(repoTitle, repoDir)
    fprintf('%s\n', repmat('=', 1, 84));
    fprintf('  %s\n', repoTitle);
    fprintf('%s\n', repmat('=', 1, 84));

    repoStats = struct('matlab', initCounter(), 'cpp', initCounter(), 'mex', initCounter());

    if ~isfolder(repoDir) || (~isfile(fullfile(repoDir, '.git')) && ~isfolder(fullfile(repoDir, '.git')))
        fprintf('Directory not found or not a git repository: %s\n\n', repoDir);
        return;
    end

    % Get git-tracked files relative to repoDir
    cmd = sprintf('git -C "%s" ls-files', repoDir);
    [status, cmdOut] = system(cmd);
    if status ~= 0
        fprintf('Failed to list git-tracked files in: %s\n\n', repoDir);
        return;
    end

    fileList = splitlines(string(cmdOut));
    fileList(fileList == "") = [];

    cExts = [".c", ".h", ".cpp", ".hpp", ".cc", ".cxx", ".hxx", ".C", ".H"];

    for i = 1:numel(fileList)
        relPath = fileList(i);
        fullPath = fullfile(repoDir, relPath);

        if ~isfile(fullPath)
            continue;
        end

        [~, ~, ext] = fileparts(fullPath);

        if strcmpi(ext, '.m')
            c = countFile(fullPath, 'matlab');
            repoStats.matlab = addCounts(repoStats.matlab, c);
        elseif any(strcmpi(ext, cExts))
            if isMexGateway(fullPath)
                c = countFile(fullPath, 'c_cpp');
                repoStats.mex = addCounts(repoStats.mex, c);
            else
                c = countFile(fullPath, 'c_cpp');
                repoStats.cpp = addCounts(repoStats.cpp, c);
            end
        end
    end

    % Compute Total across all categories
    totalRepo = addCounts(addCounts(repoStats.matlab, repoStats.cpp), repoStats.mex);

    % Print single unified table
    printUnifiedTable(repoStats, totalRepo);
end

function isMex = isMexGateway(filePath)
    isMex = false;
    fid = fopen(filePath, 'r');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid));

    content = fread(fid, '*char')';
    if ~isempty(regexp(content, '\bmexFunction\b', 'once')) || ...
       contains(content, 'matlab::mex::Function')
        isMex = true;
    end
end

function c = countFile(filePath, langType)
    c = struct('files', 1, 'comment', 0, 'code', 0);
    fid = fopen(filePath, 'rt', 'native', 'UTF-8');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid));

    inBlockComment = false;

    if strcmp(langType, 'matlab')
        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line)
                break;
            end
            str = strtrim(line);
            if isempty(str)
                continue;
            end

            if startsWith(str, '%{')
                inBlockComment = true;
                c.comment = c.comment + 1;
                continue;
            end
            if startsWith(str, '%}') && inBlockComment
                inBlockComment = false;
                c.comment = c.comment + 1;
                continue;
            end

            if inBlockComment || startsWith(str, '%')
                c.comment = c.comment + 1;
            else
                c.code = c.code + 1;
            end
        end
    else
        % C / C++
        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line)
                break;
            end
            str = strtrim(line);
            if isempty(str)
                continue;
            end

            if inBlockComment
                c.comment = c.comment + 1;
                if contains(str, '*/')
                    inBlockComment = false;
                end
                continue;
            end

            if startsWith(str, '/*')
                c.comment = c.comment + 1;
                if ~contains(str, '*/')
                    inBlockComment = true;
                end
                continue;
            end

            if startsWith(str, '//')
                c.comment = c.comment + 1;
            else
                c.code = c.code + 1;
            end
        end
    end
end

function cnt = initCounter()
    cnt = struct('files', 0, 'comment', 0, 'code', 0);
end

function total = addCounts(total, c)
    total.files   = total.files + c.files;
    total.comment = total.comment + c.comment;
    total.code    = total.code + c.code;
end

function printUnifiedTable(repoStats, totalRepo)
    fprintf('%-30s %7s %17s %17s %9s\n', 'Category', 'files', 'comment', 'code', 'total');
    fprintf('%s\n', repmat('-', 1, 84));

    printRow('MATLAB (.m)', repoStats.matlab);
    printRow('C / C++ (sources & headers)', repoStats.cpp);
    printRow('MEX Gateways (entry points)', repoStats.mex);

    fprintf('%s\n', repmat('-', 1, 84));
    printRow('TOTAL', totalRepo);
    fprintf('\n');
end

function printRow(label, cnt)
    totalLines = cnt.comment + cnt.code;
    if totalLines > 0
        pctCode = (cnt.code / totalLines) * 100;
        pctComment = (cnt.comment / totalLines) * 100;
    else
        pctCode = 0;
        pctComment = 0;
    end

    strComment = sprintf('%d (%.1f%%)', cnt.comment, pctComment);
    strCode    = sprintf('%d (%.1f%%)', cnt.code, pctCode);

    fprintf('%-30s %7d %17s %17s %9d\n', label, cnt.files, strComment, strCode, totalLines);
end
