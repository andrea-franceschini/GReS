%% benchmark_ArrayOfArrays_vs_cell.m
% Compare row-by-row access speed between:
%   1) ArrayOfArrays.getRowsMatrix(r)
%   2) Standard cell-array access C{r}
%
% The benchmark tests both sequential and randomized row access.
% A checksum is accumulated so MATLAB cannot discard the accessed data.

clear
clc

%% User parameters
nRows       = 1e5;
minRowLen   = 2;
maxRowLen   = 12;
nRepeats    = 3;
rngSeed     = 1;

%% Generate test data
rng(rngSeed)

rowLengths = randi([minRowLen,maxRowLen],nRows,1);

C = cell(nRows,1);
for r = 1:nRows
  C{r} = rand(1,rowLengths(r));
end

% Construct ArrayOfArrays from flattened data and row lengths.
flat = [C{:}];
A = ArrayOfArrays(flat,rowLengths);

%% Validate that both containers return identical rows
testRows = unique(randi(nRows,min(1000,nRows),1));

for k = 1:numel(testRows)
  r = testRows(k);

  rowA = A.getRowsMatrix(r);
  rowC = C{r};

  assert(isequal(rowA,rowC), ...
    'Mismatch between ArrayOfArrays and cell array at row %d.',r);
end

%% Access orders
sequentialOrder = (1:nRows)';
randomOrder     = randperm(nRows)';

%% Warm-up
checksumASequential(A,sequentialOrder);
checksumCSequential(C,sequentialOrder);
checksumARandom(A,randomOrder);
checksumCRandom(C,randomOrder);

%% Benchmark
tASequential = zeros(nRepeats,1);
tCSequential = zeros(nRepeats,1);
tARandom     = zeros(nRepeats,1);
tCRandom     = zeros(nRepeats,1);

for k = 1:nRepeats
  tASequential(k) = timeit(@() checksumASequential(A,sequentialOrder));
  tCSequential(k) = timeit(@() checksumCSequential(C,sequentialOrder));

  tARandom(k) = timeit(@() checksumARandom(A,randomOrder));
  tCRandom(k) = timeit(@() checksumCRandom(C,randomOrder));
end

%% Results
result = table( ...
  ["ArrayOfArrays"; "Cell array"; "ArrayOfArrays"; "Cell array"], ...
  ["Sequential"; "Sequential"; "Random"; "Random"], ...
  [median(tASequential); median(tCSequential); ...
   median(tARandom); median(tCRandom)], ...
  'VariableNames',{'Container','AccessOrder','MedianTime_s'});

result.TimePerRow_ns = 1e9 .* result.MedianTime_s ./ nRows;

disp(result)

fprintf('\nSequential access:\n')
fprintf('  ArrayOfArrays / cell time ratio: %.3f\n', ...
  median(tASequential)/median(tCSequential));

fprintf('\nRandom access:\n')
fprintf('  ArrayOfArrays / cell time ratio: %.3f\n', ...
  median(tARandom)/median(tCRandom));

fprintf('\nMemory usage:\n')
infoA = whos('A');
infoC = whos('C');
fprintf('  ArrayOfArrays object variable: %.3f MB\n',infoA.bytes/2^20);
fprintf('  Cell array variable:           %.3f MB\n',infoC.bytes/2^20);
fprintf(['  Note: whos may not report all memory referenced internally by ', ...
         'a handle object.\n']);

%% Optional plot
figure
bar(categorical( ...
  result.Container + " - " + result.AccessOrder), ...
  result.TimePerRow_ns)

ylabel('Time per row [ns]')
title(sprintf('Row access benchmark, %d rows',nRows))
grid on

%% Local benchmark functions

function checksum = checksumASequential(A,order)
  checksum = 0;

  for k = 1:numel(order)
    row = A.getRowsMatrix(order(k));

    % Use multiple values so the entire row is materially accessed.
    checksum = checksum + sum(row);
  end
end

function checksum = checksumCSequential(C,order)
  checksum = 0;

  for k = 1:numel(order)
    row = C{order(k)};
    checksum = checksum + sum(row);
  end
end

function checksum = checksumARandom(A,order)
  checksum = 0;

  for k = 1:numel(order)
    row = A.getRowsMatrix(order(k));
    checksum = checksum + sum(row);
  end
end

function checksum = checksumCRandom(C,order)
  checksum = 0;

  for k = 1:numel(order)
    row = C{order(k)};
    checksum = checksum + sum(row);
  end
end
