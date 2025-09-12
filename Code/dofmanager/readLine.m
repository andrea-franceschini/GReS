function token = readLine(fID, fName)
  if feof(fID)
    error('Unexpected end of file %s before End statement',fName);
  end
  token = sscanf(fgetl(fID), '%s');
end