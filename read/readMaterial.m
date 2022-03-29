function [numMat,MatData] = readMaterial(fileName)
%function for materials' input data

%opening input file:
fid = fopen(fileName,'r');

%reading total number of materials as a cell:
numMat = textscan(fid,'%f %*s',1);
%converting cell array into ordinary array:
numMat = cell2mat(numMat);

%reading each materials' data in a block
%initializing the first material:
Block = 1;

%reading until the end of the input file:
while (~feof(fid))
%For each material (Block):

%reading material's ID:
MatID = textscan(fid,'%f %*s',1);
%celldisp(MatID);

%reading material's mathematical model as a TAG:
MatModel = textscan(fid,'%f %*s',1);
%celldisp(MatModel);

%reading material's data in double precision, specifying that the data after
%the "!" sign is to be ignored
InputData = textscan(fid,'%f','CommentStyle','!');
%converting cell array into ordinary array:
InputData = cell2mat(InputData);
%defining the number of columns of output matrix "MatData":
nCols = size(InputData);
nCols = nCols(1)+2;
%converting InputData into a row vector:
InputData = InputData';
%disp(InputData)

%read and discard end of the Block:
eof = textscan(fid,'%s',1,'Delimiter','\n');

% ------ Creating matrix "MatData" for materials' mechanical data ------

%assigning MatID in the first column:
MatData(Block,1) = cell2mat(MatID);

%assigning MatModel in the second column:
MatData(Block,2) = cell2mat(MatModel);

%assigning all material's properties in the remaining columns:
MatData(Block,3:nCols) = InputData;
%disp(MatData)

%increase Block index:
Block = Block+1;

end

end

