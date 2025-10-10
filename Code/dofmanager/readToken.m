% Check for eof and read the next token
function token = readToken(fID, fName)
  if feof(fID)
    error('Unexpected end of file %s before End statement',fName);
  end
  token = sscanf(fgetl(fID), '%s',1);
end