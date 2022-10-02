classdef Materials < handle
  % MATERIAL General material class

  properties (Access = public)
    % Creation of a Map object 
    db = containers.Map;
  end

  methods (Access = public)
    % Class constructor method   
    function obj = Materials(fileName)
      % Calling the function to read input data from file
      obj.readInputFile(fileName)
    end

    % Check if the matIdentifier defined by the user is a key of the Map object
    function mat = getMaterial(obj,matIdentifier)
      if (obj.db.isKey(matIdentifier))
        mat = obj.db(matIdentifier);
      else
        % Displaying error message if the material class has not been created yet
        error('Material % not present', matIdentifier);
      end
    end
  end

  methods (Access = private)
    % Reading material input file
    function readInputFile(obj, fileName)
      fid = fopen(fileName, 'r');
      while (~feof(fid))
        line = obj.getNewLine(fid);
        if (length(line) == 0)
          continue;
        end
        % Reading material name (first line)
        matName = sscanf(line,'%s',1);
        line = obj.getNewLine(fid);
        % Reading material identifier (second line)
        matIdentifier = sscanf(line,'%s',1);
        data = [];
        % Reading material data until char = End
        while (~strcmp(line, 'End'))
          line = obj.getNewLine(fid);
          parts = strsplit(line, '%');
          line = parts{1};
          if (~strcmp(line, 'End'))
            data = [data,line];
          end
        end
        
       % Calling the specific material class based on matName
        switch lower(matName)
          case 'elastic'
            obj.db(matIdentifier) = Elastic(data);
          case 'hypoelastic'
            obj.db(matIdentifier) = HypoElastic(data);
          case 'transvelastic'
            obj.db(matIdentifier) = TransvElastic(data);
          case 'porousrock'
            obj.db(matIdentifier) = PorousRock(data);
          case 'fluid'
            obj.db(matIdentifier) = Fluid(data);
          otherwise
            error('%s not available', matName);
        end
      end
      fclose(fid);
    end

    % Reading lines from the input file until the end of the file
    function line = getNewLine(obj, fid)
      line = fgetl(fid);
      while (~feof(fid) && length(line) == 0)
        line = fgetl(fid);
      end
    end

  end
end
