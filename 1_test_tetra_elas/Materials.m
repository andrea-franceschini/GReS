classdef Materials < handle
  % MATERIAL General material class

  properties (Access = public)
    % Creation of a Map Object 
    db = containers.Map;
  end


  methods (Access = public)
    % Class constructor method   
    function obj = Materials(fileName)
      % Calling the function readInputFile to read input data from file
      obj.readInputFile(fileName)
    end

    % Function to get the matIdentifier defined by the user 
    function mat = getMaterial(obj,matIdentifier)
     % Determining if the matIdentifier is related to one of the existing 
     % materials' class
      if (obj.db.isKey(matIdentifier))
        mat = obj.db(matIdentifier);
      else
        % Displaying error message if the material class has not been
        % created yet
        error('Material % not present', matIdentifier);
      end
    end
  end

  methods (Access = private)
      % Function to read the material input file:
    function readInputFile(obj, fileName)
      fid = fopen(fileName, 'r');
      % Reading until the end of the file:
      while (~feof(fid))
        % Reading the first line of the input file
        line = obj.getNewLine(fid);
        % Reading each line of the file until it's empty
        if (length(line) == 0)
          continue;
        end
        % Reading material's name (first line)
        matName = sscanf(line,'%s',1);
        % Reading the second line of the input file
        line = obj.getNewLine(fid);
        % Reading material's identifier
        matIdentifier = sscanf(line,'%s',1);
        % Creation of an empty array to store data
        data = [];
        % Reading material's data until line = End
        while (~strcmp(line, 'End'))
          line = obj.getNewLine(fid);
          parts = strsplit(line, '%');
          line = parts{1};
          if (~strcmp(line, 'End'))
             % Filling the empty array with material's data
            data = [data,line];
          end
        end
        
       % Calling the right material class based on the matName
       % Using lower to convert "matName" characters in lowercase character
        switch lower(matName)
          case 'elastic'
            obj.db(matIdentifier) = Elastic(data);
          case 'hypoplastic'
            obj.db(matIdentifier) = HypoPlastic(data);
          case 'hypoelastic'
            obj.db(matIdentifier) = HypoElastic(data);
          case 'transvelastic'
            obj.db(matIdentifier) = TransvElastic(data);
%           case 'camclay'
%             obj.db(matIdentifier) = ModCamClay(data);
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

    % Function for reading lines from the input file until the end of the
    % file
    function line = getNewLine(obj, fid)
      line = fgetl(fid);
      while (~feof(fid) && length(line) == 0)
        line = fgetl(fid);
      end
    end

  end
end
