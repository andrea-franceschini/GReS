classdef Materials < handle
  % MATERIAL - General material class

  properties (Access = private)
    % Creation of a Map object 
    db = containers.Map('KeyType','double','ValueType','any');
  end

  methods (Access = public)
    % Class constructor method   
    function obj = Materials(fileName)
      % Calling the function to read input data from file
      obj.readInputFile(fileName)
    end

    % Get the material defined by matIdentifier and check if it is a
    % key of the Map db
    function mat = getMaterial(obj,matID)
      %
      % The preliminary check whether matID is key of db has been commented
      % since it is highly expensive
%       if (obj.db.isKey(matID))
        mat = obj.db(matID);
%       else
%       Displaying error message if matIdentifier is not a key 
%       of the map
%         error('Material %s not present', matID);
%       end
    end
    
    %
    % Destructor
    function delete(obj)
      remove(obj.db,keys(obj.db));
    end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, fileName)
      fid = fopen(fileName, 'r');
      % flBlock: flag reporting the reading status
      %          0 -> the material block has been read correctly
      %          1 -> error while reading the material block
      %               (end-of-file before the 'End' statement
      %          2 -> end of file reached
      flBlock = 0;
      % Number of blocks counter
      nBlock = 0;
      while (flBlock == 0)
        % Update the counter
        nBlock = nBlock + 1;
        % Reading the material block
        [flBlock,block] = Materials.readBlock(fid);
        if flBlock == 0
          if ~isnan(str2double(block(1)))
              error('The first entry in block %d of Materials file %s must be strings',nBlock,fileName);
          end
          % Calling the specific material class based on the
          % material name
          switch lower(block(1))
            case 'elastic'
              obj.db(nBlock) = Elastic(block);
            case 'hypoelastic'
              obj.db(nBlock) = HypoElastic(block);
            case 'transvelastic'
              obj.db(nBlock) = TransvElastic(block);
            case 'porousrock'
              obj.db(nBlock) = PorousRock(block);
            case 'fluid'
              obj.db(nBlock) = Fluid(block);
            otherwise
              error('Material %s not available (block %d of Materials)', block(1),nBlock);
          end
        elseif flBlock == 1
          error('Error encountered while reading block %d of Material file',nBlock);
        end
      end
      fclose(fid);
    end
  end

  methods (Static = true)
    % Read the material block
    function [flBlock,block] = readBlock(fid)
      block = [];
      % flBegBlock: flag marking the beginning of a material block
      flBegBlock = false;
      [flEof,line] = Materials.readLine(fid);
      % Dealing with the initial blank lines, if any
      % flEof: end-of-file flag
      %        1 -> end of file reached
      %        0 -> otherwise
      while (isempty(strtrim(line)) && flEof == 0)
        [flEof,line] = Materials.readLine(fid);
      end
      %
      if flEof == 0
        flBegBlock = true;
      end
      % Reading the material block
      while (~strcmp(line,'End') && flEof == 0)
        line = string(strtrim(strtok(line,'%')));
        block = [block; line];
        [flEof,line] = Materials.readLine(fid);
      end
      %
      if flEof == 0 && flBegBlock
        flBlock = 0;   % block read correctly
      elseif flEof == 1 && flBegBlock
        flBlock = 1;   % end of file while reading
      elseif flEof == 1 && ~flBegBlock
        flBlock = 2;   % end of file reached
      end
    end

    % Read the next line and check for eof
    function [flEof,line] = readLine(fid)
      flEof = feof(fid);   % end-of-file flag
      if flEof == 1
        line = '';
      else
        line = strtrim(fgetl(fid));
      end  
    end
  end
end