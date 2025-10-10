function [nVals, vals] = readEntitySet(fileName)
  if (~exist(fileName, 'file'))
    error('File %s does not seem to exist. Please, check the provided file.', fileName);
  end
  header = false;
  fid = fopen(fileName, 'r');
  while (~feof(fid) && ~header)
    line = fgetl(fid);
    word = sscanf(line, '%s');
    if (~strcmp(word(1), '%'))
      % If this is not a commented line (not starting with %)
      nVals = sscanf(line, '%i');
      header = true;
    end
  end
  if (~header)
    error('Missing header in readSet.');
  end
  nValMax = sum(nVals);
  vals = zeros(nValMax,1);
  id = 1;
  while ~feof(fid)
    line = fgetl(fid);
    word = sscanf(line, '%s');
    if (~strcmp(word(1), '%'))
      % If this is not a commented line (not starting with %)
      num = sscanf(line, '%i');
      nNum = length(num);
      vals(id:id+nNum-1) = num;
      id = id + nNum;
    end
  end
  fclose(fid);
end
