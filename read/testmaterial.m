% ------- Creation of material INPUT file "MaterialINPUT.mat" --------

% fid = fopen('MatINPUT.mat','w');
% fprint(fid,'%s \n', '%MaterialNames')
% fprint(fid,'%d \t %s \t %s \n',ID,MatModel,Name)
% fprint(fid,'%s \n %s \n %s \n','%EndMaterialNames','%Data','%Mat 1')
% fprint(fid,'%s \n', '%Data')

close all;
clear;
clc;

%fileName = 'matINPUT.mat';
fileName = uigetfile('*.mat');

% %READMATERIAL for reading materials' input data
% 
% %opening input file
% fid = fopen(fileName,'r');
% 
% %reading total number of materials as a cell
% numMat = textscan(fid,'%f %*s',1);
% numMat = cell2mat(numMat);
% 
% %reading each materials' data in a block
% %initializing first material
% Block = 1;
% 
% %reading until the end of the input file
% while (~feof(fid))
% %For each material (Block):
% 
% %reading material's ID
% MatID = textscan(fid,'%f %*s',1);
% %celldisp(MatID);
% 
% %reading material's mathematical model as a TAG
% MatModel = textscan(fid,'%f %*s',1);
% %celldisp(MatModel);
% 
% %reading material's data
% InputData = textscan(fid,'%f','CommentStyle','!');
% InputData = cell2mat(InputData);
% nCols = size(InputData);
% nCols = nCols(1)+2;
% InputData = InputData';
% disp(InputData)
% 
% %read and discard end of the Block
% eof = textscan(fid,'%s',1,'Delimiter','\n');
% 
% %Creating matrix "MatData" for materials' general info:
% 
% %assigning MatID in the first column
% MatData(Block,1) = cell2mat(MatID);
% 
% %assigning MatModel in the second column
% MatData(Block,2) = cell2mat(MatModel);
% 
% %assigning material's properties in the remaining column
% MatData(Block,3:nCols) = InputData;
% disp(MatData)
% 
% %increase Block index
% Block = Block+1;
% 
% end

mat = Material ();
tic;
[numMat,MatData] = readMaterial(fileName);
setMaterial(mat,fileName)
t1 = toc;
fprintf('Time to read %.3f [s]\n', t1)


