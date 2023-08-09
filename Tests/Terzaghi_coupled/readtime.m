function [t] = readtime(fileName)
fid = fopen(fileName,'r');

%get tag for end of file reached
flEof = feof(fid);
if flEof ==1
    line = '';
else
    line = strtrim(fgetl(fid)); 
end
block = '';

while ~strcmp(line,'End')
    line = strtrim(line);
    if isempty(line)
        error('blank line encountered')
    elseif flEof == 1 
        error('end of file reached')
    end
    block = [block, line, ' '];
    flEof = feof(fid);
    if flEof ==1
        line = '';
    else
        line = strtrim(fgetl(fid)); 
    end
end

blockSplt = strsplit(string(deblank(block)));
t = str2double(blockSplt);
    
        


