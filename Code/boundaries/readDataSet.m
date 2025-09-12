function vals = readDataSet(fileName, nVals)
  if (~exist(fileName, 'file'))
    error('File %s does not seem to exist. Please, check the provided file.', fileName);
  end
  fid = fopen(fileName, 'r');
  vals = zeros(nVals,1);
  id = 1;
  while ~feof(fid)
    line = fgetl(fid);
    word = sscanf(line, '%s');
    if (~strcmp(word(1), '%'))
      % If this is not a commented line (not starting with %)
      num = sscanf(line, '%e');
      nNum = length(num);
      vals(id:id+nNum-1) = num;
      id = id + nNum;
    end
  end
  if length(vals) ~= nVals
      error('Number of values in %s not matching number of constrained entities',fileName)
  end
  fclose(fid);
end
