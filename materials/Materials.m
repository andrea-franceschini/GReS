classdef Materials < handle
  %MATERIAL General material class

  properties (Access = private)
    db = containers.Map;
  end

  methods (Access = public)
    function obj = Materials(fileName)
      obj.readInputFile(fileName);
    end

    function mat = getMaterial(obj, matIdentifier)
      if (obj.db.isKey(matIdentifier))
        mat = obj.db(matIdentifier);
      else
        error('Material % not present', matIdentifier);
      end
    end
  end

  methods (Access = private)
    function readInputFile(obj, fileName)
      fid = fopen(fileName, 'r');
      while (~feof(fid))
        line = obj.getNewLine(fid);
        if (length(line) == 0)
          continue;
        end
        matName = sscanf(line,'%s',1);
        line = obj.getNewLine(fid);
        matIdentifier = sscanf(line,'%s',1);
        data = [];
        while (~strcmp(line, 'End'))
          line = obj.getNewLine(fid);
          parts = strsplit(line, '%');
          line = parts{1};
          if (~strcmp(line, 'End'))
            data = [data,line];
          end
        end
        switch lower(matName)
          case 'elastic'
            obj.db(matIdentifier) = Elastic(data);
          case 'hypoplastic'
            obj.db(matIdentifier) = HypoPlastic(data);
          case 'porousrock'
            obj.db(matIdentifier) = PorousRock(data);
          otherwise
            error('%s not available', matName);
        end
      end
      fclose(fid);
    end

    function line = getNewLine(obj, fid)
      line = fgetl(fid);
      while (~feof(fid) && length(line) == 0)
        line = fgetl(fid);
      end
    end

  end
end
