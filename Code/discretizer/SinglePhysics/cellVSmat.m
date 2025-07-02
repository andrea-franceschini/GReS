% Number of "keys"
N = 3e4;  % 1 million entries

% Randomly generate array sizes (between 1 and 4 values per key)
valueCounts = randi(4, 1, N);

% -------------------------
% 1. CELL ARRAY
% -------------------------
cellMap = cell(1, N);
for i = 1:N
    cellMap{i} = randi(100, 1, valueCounts(i));
end

% -------------------------
% 2. SPARSE MATRIX
% -------------------------
% We'll store all values in a sparse matrix with up to 4 columns
sparseMap = sparse(N, 4);
for i = 1:N
    sparseMap(i, 1:valueCounts(i)) = randi(100, 1, valueCounts(i));
end

% -------------------------
% 3. STANDARD ARRAY
% -------------------------
% Flatten all values into a single array with zeros between
maxEntries = sum(valueCounts);
flatArray = zeros(1, maxEntries);
offsets = zeros(1, N+1)+1;
current = 1;

for i = 1:N
    nvals = valueCounts(i);
    flatArray(current:current+nvals-1) = randi(100, 1, nvals);
    offsets(i+1) = current + nvals;
    current = current + nvals;
end

% -------------------------
% ACCESS BENCHMARK
% -------------------------

% 1. Access cell array
tic
total1 = 0;
for i = 1:N
    data = cellMap{i};
    total1 = total1 + sum(data);  % simulate access
end
time_cell = toc;
disp(['Cell array access time: ', num2str(time_cell), ' s']);

% 2. Access sparse matrix
tic
total2 = 0;
for i = 1:N
    data = nonzeros(sparseMap(i, :))';
    total2 = total2 + sum(data);
end
time_sparse = toc;
disp(['Sparse matrix access time: ', num2str(time_sparse), ' s']);

% 3. Access flattened array using offsets
tic
total3 = 0;
for i = 1:N
    idx1 = offsets(i);
    idx2 = offsets(i+1) - 1;
    data = flatArray(idx1:idx2);
    total3 = total3 + sum(data);
end
time_flat = toc;
disp(['Flat array access time: ', num2str(time_flat), ' s']);

% -------------------------
% Compare Results
% -------------------------
fprintf('\nSummary:\n');
fprintf('Cell   : %.3f s\n', time_cell);
fprintf('Sparse : %.3f s\n', time_sparse);
fprintf('Flat   : %.3f s\n', time_flat);

speed_cell = time_cell / time_flat;
speed_sparse = time_sparse / time_flat;

fprintf('\nFlat array is:\n');
fprintf(' - %.2fx faster than cell array\n', speed_cell);
fprintf(' - %.2fx faster than sparse matrix\n', speed_sparse);
