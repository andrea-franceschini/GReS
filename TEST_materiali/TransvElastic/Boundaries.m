classdef Boundaries < handle
  % BOUNDARY CONDITIONS General boundary conditions class

  properties (Access = private)
    % Creation of a Map object 
    db = containers.Map;
  end

  methods (Access = public)
    % Class constructor method
    function obj = Boundaries(fileName)
      % Calling the function to read input data from file
      obj.readInputFile(fileName)
    end

    % Check if the BCIdentifier defined by the user is a key of the Map object
    function bcType = getBC(obj, BCIdentifier)
      if (obj.db.isKey(BCIdentifier))
        bcType = obj.db(BCIdentifier);
      else
        % Displaying error message if the boundary class has not been created yet
         error('Boundary condition % not present', BCIdentifier);
      end
    end
  end

  methods (Access = private)
    % Reading boundary input file
    function readInputFile(obj, fileName)
      fid = fopen(fileName, 'r');
      while (~feof(fid))
        line = obj.getNewLine(fid);
        if (length(line) == 0)
          continue;
        end
        % Reading BC name (first line)
        BCName = sscanf(line,'%s',1);
        line = obj.getNewLine(fid);
        % Reading BC identifier (second line)
        BCIdentifier = sscanf(line,'%s',1);
        data = [];
        % Reading BC data until char = End 
        while (~strcmp(line, 'End'))
          line = obj.getNewLine(fid);
          parts = strsplit(line,'%');
          line = parts{1};
          if (~strcmp(line, 'End'))
            data = [data,line];
          end
        end        
       
        % Calling the specific boundary condition class based on  BCName
        switch lower(BCName)
          case 'nodebc'
            obj.db(BCIdentifier) = NodeBC(data);
          otherwise
            error('%s not available', BCName);
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









