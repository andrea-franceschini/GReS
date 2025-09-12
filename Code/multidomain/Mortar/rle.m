function [outMat,count] = rle(inMat)
% the input matrix has already been sorted by rows
% Find unique rows and their first occurrences
[outMat, firstIdx, ~] = unique(inMat, 'rows', 'stable');

% Compute counts of each unique row
count = diff([firstIdx; size(inMat,1) + 1]);

% Restore original order of indices
%indices = sortIdx(firstIdx);
end