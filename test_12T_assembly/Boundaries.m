classdef Boundaries < handle
  % BOUNDARIES General material class

  properties (Access = private)
    % Creation of a Map Object 
    db = containers.Map;
  end

  methods (Access = public)
      % Class constructor method
    function obj = Boundaries(fileName)
      % Reading input data from file
      obj.readInputFile(fileName)
    end

    % Getting BCidentifier defined by the user
    function bcType = getBC(obj, BCIdentifier)
        % Determining if the BCIdentifier is related to one of the
        % existing boundaries' class
      if (obj.db.isKey(BCIdentifier))
        bcType = obj.db(BCIdentifier);
      else
        %displaying error message if the boundary class has not been
        %created yet
         error('Boundary condition % not present', BCIdentifier);
      end
    end
  end

  methods (Access = private)
      
    %Function to read the input file
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
        % Reading BC's name (first line)
        BCName = sscanf(line,'%s',1);
        % Reading the second line of the input file
        line = obj.getNewLine(fid);
        % Reading boundary condition's identifier
        BCIdentifier = sscanf(line,'%s',1);
        % Creation of an empty array to store data
        data = [];
        % Reading BC's data until line = End
        while (~strcmp(line, 'End'))
          line = obj.getNewLine(fid);
          parts = strsplit(line,'%');
          line = parts{1};
          if (~strcmp(line, 'End'))
             % Filling the empty array with material's data
            data = [data,line];
          end
        end        
       
        
  % Calling the right BC class based on BCIdentifier 
        switch lower(BCName)
          case 'nodaldirichlet'
            obj.db(BCIdentifier) = NodeBC(data);
          case 'nodalneumann'
            obj.db(BCIdentifier) = NodeBC(data);
%           case 'superficialdirichlet'
%             obj.db(BCIdentifier) = FaceBC(data);
%           case 'superficialneumann'
%             obj.db(BCIdentifier) = FaceBC(data);
          otherwise
            error('%s not available', matName);
        end
       end
      fclose(fid);
    end


    
    % Function to read lines from the input file until the end of the file
    function line = getNewLine(obj, fid)
      line = fgetl(fid);
      while (~feof(fid) && length(line) == 0)
        line = fgetl(fid);
      end
    end

  end
end









