% Check for eof and read the next token
function [tags,token] = readTokenList(fID, fName)
  if feof(fID)
    error('Unexpected end of file %s before End statement',fName);
  end
  line = strtok(fgetl(fID),'%');
  [tags,~,~,n] = sscanf(line, '%i');
  token = sscanf(line(n:end), '%s', 1);
end