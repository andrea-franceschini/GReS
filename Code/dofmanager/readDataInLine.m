function vec = readDataInLine(fID, matFileName, nEntry)
  if feof(fID)
    error('Unexpected end of file %s while reading',matFileName);
  else
    [vec, n, errmsg] = sscanf(fgetl(fID), '%e', nEntry);
    if n ~= nEntry
      error('Wrong number of row entries in material file %s',matFileName);
    elseif ~isempty(errmsg)
      error('%s in file %s', errmsg, matFileName);
    end
  end
end