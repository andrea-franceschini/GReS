function compare_solver_stats(varargin)
    % COMPARE_SOLVER_STATS Plots a 2x2 performance dashboard for N solver log files
    % and prints a statistical summary table for setup times, iterations, and total runtimes.
    % 
    % Usage:
    %   compare_solver_stats('file1.dat', 'Label 1', 'file2.dat', 'Label 2', ...)
    %   compare_solver_stats('file1.dat', 'Label 1', ..., 'saveFigs', true)
    
    if nargin < 1
        error('At least one file path must be provided.');
    end

    % --- 0. Parse 'saveFigs' Parameter ---
    saveFigs = false;
    cleanVarargin = {};
    i = 1;
    while i <= length(varargin)
        if (ischar(varargin{i}) || isstring(varargin{i})) && strcmpi(varargin{i}, 'saveFigs')
            if i + 1 <= length(varargin)
                saveFigs = logical(varargin{i+1});
                i = i + 2;
                continue;
            else
                error('The "saveFigs" option must be followed by a boolean (true or false).');
            end
        end
        cleanVarargin{end+1} = varargin{i}; %#ok<AGROW>
        i = i + 1;
    end

    if isempty(cleanVarargin)
        error('At least one file path must be provided.');
    end

    % --- 1. Parse Input Pairs (fileName, label) ---
    files = {};
    labels = {};
    
    idx = 1;
    fileCount = 1;
    while idx <= length(cleanVarargin)
        files{end+1} = char(cleanVarargin{idx}); %#ok<AGROW>
        
        % Check if a label string was explicitly provided for this file
        if idx + 1 <= length(cleanVarargin) && (ischar(cleanVarargin{idx+1}) || isstring(cleanVarargin{idx+1}))
            labels{end+1} = char(cleanVarargin{idx+1}); %#ok<AGROW>
            idx = idx + 2;
        else
            labels{end+1} = sprintf('Dataset %d', fileCount); %#ok<AGROW>
            idx = idx + 1;
        end
        fileCount = fileCount + 1;
    end
    numFiles = length(files);

    % --- 2. Load Data for All Files & Pool Times ---
    data = cell(1, numFiles);
    allPrecTimes = [];
    for k = 1:numFiles
       if ~exist(files{k}, 'file')
          error('Unable to find or open file: %s', files{k});
       end
       T = readSolverLog(files{k});
       s.x = T.SolN;
       s.iter = T.Iter;
       s.solTime = T.SolTime;
       s.precTime = T.PrecTime;
       s.cumTime = cumsum(s.precTime + s.solTime);
       data{k} = s;
       % Accumulate non-zero setup times across all files
       allPrecTimes = [allPrecTimes; s.precTime(s.precTime > 0)]; %#ok<AGROW>
    end

    % --- Automatically Compute Separation Threshold ---
    timeThreshold = findSetupThreshold(allPrecTimes);
    fprintf('Auto-detected SAM/AMG setup time cutoff: %.4f s\n', timeThreshold);

    % --- Compute Summary Statistics per Dataset ---
    nAMG = zeros(numFiles, 1);
    avgTimeAMG = zeros(numFiles, 1);
    nSAM = zeros(numFiles, 1);
    avgTimeSAM = zeros(numFiles, 1);
    avgIter = zeros(numFiles, 1);
    avgSolTime = zeros(numFiles, 1);
    totalTime = zeros(numFiles, 1);

    for k = 1:numFiles
       s = data{k};
       maskAMG = s.precTime >= timeThreshold;
       maskSAM = s.precTime > 0 & s.precTime < timeThreshold;
       nAMG(k) = sum(maskAMG);
       if nAMG(k) > 0
          avgTimeAMG(k) = mean(s.precTime(maskAMG));
       else
          avgTimeAMG(k) = NaN;
       end
       nSAM(k) = sum(maskSAM);
       if nSAM(k) > 0
          avgTimeSAM(k) = mean(s.precTime(maskSAM));
       else
          avgTimeSAM(k) = NaN;
       end
       % Overall Iteration, Solve Time Averages, & Total Time
       avgIter(k) = mean(s.iter);
       avgSolTime(k) = mean(s.solTime);
       totalTime(k) = sum(s.precTime) + sum(s.solTime);
    end

    % --- 3. Print Summary Table to Command Window ---
    summaryTable = table(labels', nAMG, avgTimeAMG, nSAM, avgTimeSAM, avgIter, avgSolTime, totalTime, ...
        'VariableNames', {'Label', 'nAMG', 'AvgTimeAMG', 'nSAM', 'AvgTimeSAM', 'AvgIter', 'AvgSolTime_s', 'TotalTime_s'});
    
    fprintf('\n=========================================== SOLVER SUMMARY METRICS =======================================\n');
    disp(summaryTable);
    fprintf('==========================================================================================================\n\n');

    % --- 4. Graphics Setup ---
    colors = lines(numFiles);
    markers = {'o', 's', '^', 'v', 'd', 'p', 'h', 'x', '+'};
    fig = figure('Color', 'w', 'Position', [100, 100, 1100, 750]);
    if numFiles > 1
        sgtitle('Performance Analysis: Solver Comparison', 'FontSize', 14, 'FontWeight', 'bold');
    else
        sgtitle(['Performance Analysis: ', labels{1}], 'FontSize', 14, 'FontWeight', 'bold');
    end
    ln = 1.5;

    % --- Subplot 1: GMRES Iterations (Semilog Y) ---
    ax1 = subplot(2, 2, 1); hold on; grid on; grid minor; box on;
    set(gca, 'YScale', 'log');
    for k = 1:numFiles
        m = markers{mod(k-1, length(markers)) + 1};
        semilogy(data{k}.x, data{k}.iter, ['-', m], 'Color', colors(k, :), 'LineWidth', ln, ...
            'MarkerSize', 1, 'MarkerFaceColor', colors(k, :), 'DisplayName', labels{k});
    end
    ylabel('Iterations (Log Scale)', 'FontSize', 16); 
    xlabel('Sequence Step', 'FontSize', 16);
    title('GMRES Iterations', 'FontSize', 18); 
    legend('Location', 'northwest', 'FontSize', 16);

    % --- Subplot 2: Setup Time (Points Only, No Lines) ---
    ax2 = subplot(2, 2, 2); hold on; grid on; grid minor; box on;
    set(gca, 'YScale', 'log');
    hPlots = gobjects(numFiles, 1);
    for k = 1:numFiles
       active = data{k}.precTime > 0;
       xActive = data{k}.x(active);
       yActive = data{k}.precTime(active);
       hPlots(k) = plot(xActive, yActive, 'o', ...
          'LineStyle', 'none', ...
          'MarkerSize', 5, ...
          'MarkerFaceColor', colors(k, :), ...
          'MarkerEdgeColor', 'none', ...
          'DisplayName', labels{k});
    end
    ylabel('Time (s, Log Scale)', 'FontSize', 16); 
    xlabel('Sequence Step', 'FontSize', 16);
    title('Setup Time (Preconditioner / SAM)', 'FontSize', 16); 
    legend(hPlots, labels, 'Location', 'best', 'FontSize', 13);

    % --- Subplot 3: GMRES Solve Time (Semilog Y) ---
    ax3 = subplot(2, 2, 3); hold on; grid on; grid minor; box on;
    set(gca, 'YScale', 'log');
    for k = 1:numFiles
        m = markers{mod(k-1, length(markers)) + 1};
        semilogy(data{k}.x, data{k}.solTime, ['-', m], 'Color', colors(k, :), 'LineWidth', ln, ...
            'MarkerSize', 1, 'MarkerFaceColor', colors(k, :), 'DisplayName', labels{k});
    end
    ylabel('Time (s, Log Scale)', 'FontSize', 16); 
    xlabel('Sequence Step', 'FontSize', 16);
    title('GMRES Solve Time', 'FontSize', 18); 
    legend('Location', 'northwest', 'FontSize', 16);

    % --- Subplot 4: Cumulative Runtime (Setup + Solve) ---
    ax4 = subplot(2, 2, 4); hold on; grid on; box on;
    for k = 1:numFiles
        m = markers{mod(k-1, length(markers)) + 1};
        plot(data{k}.x, data{k}.cumTime, ['-', m], 'Color', colors(k, :), 'LineWidth', ln, ...
            'MarkerSize', 1, 'MarkerFaceColor', colors(k, :), 'DisplayName', labels{k});
    end
    ylabel('Total Accumulated Time (s)', 'FontSize', 16); 
    xlabel('Sequence Step', 'FontSize', 16);
    title('Cumulative Time (Setup + Solve)', 'FontSize', 18); 
    legend('Location', 'northwest', 'FontSize', 16);

    % --- 5. Export to PDF if Requested ---
    if saveFigs
        axesHandles = {ax1, ax2, ax3, ax4};
        fileNames = { ...
            'gmres_iterations.pdf', ...
            'setup_time.pdf', ...
            'gmres_solve_time.pdf', ...
            'cumulative_time.pdf' ...
        };

        for j = 1:4
            exportSubplotAsPdf(axesHandles{j}, fileNames{j});
            fprintf('Saved: %s\n', fullfile(pwd, fileNames{j}));
        end

        % Save the full 2x2 dashboard as well
        exportgraphics(fig, 'solver_comparison_dashboard.pdf', 'ContentType', 'vector');
        fprintf('Saved: %s\n', fullfile(pwd, 'solver_comparison_dashboard.pdf'));
    end
end

% --- Helper Function: Exports a Subplot and its Legend as a Clean Standalone PDF ---
function exportSubplotAsPdf(sourceAx, outputPdf)
hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Check if the source axes has an active legend
leg = sourceAx.Legend;

if ~isempty(leg) && isvalid(leg)
   % Copy BOTH the axes and the legend simultaneously
   % Passing both handles to copyobj preserves their link
   copiedObjects = copyobj([sourceAx, leg], hFig);
   newAx = copiedObjects(1);
else
   % Axes has no legend; copy axes only
   newAx = copyobj(sourceAx, hFig);
end

% Adjust axes bounds so labels, titles, and legends fit comfortably
set(newAx, 'Units', 'normalized', 'Position', [0.15, 0.15, 0.75, 0.75]);

% Export vector graphic to PDF
exportgraphics(hFig, outputPdf, 'ContentType', 'vector');
close(hFig);
end

% --- Helper Function: Parses Pipe-Delimited Log File ---
function T = readSolverLog(filePath)
    lines = readlines(filePath);
    
    isDataOrHeader = ~startsWith(strtrim(lines), '-') & strtrim(lines) ~= "";
    cleanLines = lines(isDataOrHeader);
    
    rawHeader = split(cleanLines(1), '|');
    rawHeader = strtrim(rawHeader(2:end-1));
    
    varNames = cellfun(@(s) strrep(strtrim(s), ' ', ''), rawHeader, 'UniformOutput', false);
    varNames = strrep(varNames, '.', ''); 
    
    dataLines = cleanLines(2:end);
    numRows = length(dataLines);
    numCols = length(varNames);
    dataMat = zeros(numRows, numCols);
    for i = 1:numRows
        rowCols = split(dataLines(i), '|');
        rowCols = strtrim(rowCols(2:end-1));
        dataMat(i, :) = str2double(rowCols)';
    end
    T = array2table(dataMat, 'VariableNames', varNames);
end

% --- Helper Function: Automatically Finds Cutoff Between Fast (SAM) and Slow (AMG) Setups ---
function thresh = findSetupThreshold(times)
    posTimes = sort(unique(times(times > 0)));
    if length(posTimes) < 2
       thresh = 0.5;
       return;
    end
    logT = log10(posTimes);
    [maxGap, maxIdx] = max(diff(logT));
    if maxGap < 0.3
       thresh = 0.5;
       return;
    end
    thresh = 10^((logT(maxIdx) + logT(maxIdx + 1)) / 2);
end