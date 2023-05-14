classdef Materials < handle
  % MATERIAL - General material class

  properties (Access = private)
    % Creation of a Map object 
    db = containers.Map('KeyType','double','ValueType','any');
  end

  methods (Access = public)
    % Class constructor method   
    function obj = Materials(fListName)
      % Calling the function to read input data from file
      obj.readInputFiles(fListName)
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
    function readInputFile(obj, matFileName, matID)
      fID = fopen(matFileName, 'r');
      Materials.checkFileOpening(fID, matFileName);
      %
      assert(~feof(fID),'No material properties have been assigned in %s',matFileName);
      propName = Materials.readToken(fID, matFileName);
      while true
        % Calling the specific material class based on the
          % material name
          switch propName
            case 'Elastic'
              matProp.ConstLaw = Elastic(fID, matFileName);
            case 'HypoElastic'
              matProp.ConstLaw = HypoElastic(fID, matFileName);
            case 'TransvElastic'
              matProp.ConstLaw = TransvElastic(fID, matFileName);
            case 'PorousRock'
              matProp.PorousRock = PorousRock(fID, matFileName);
            case 'CapillaryCurve'
              matProp.CapillaryCurve = CapillaryCurve(fID, matFileName);
            case 'Fluid'
              matProp = Fluid(fID, matFileName);
            otherwise
              error('Material property %s not available in file %s', token,matFileName);
          end
%           [str, flEoF] = Materials.readToken(fID);
          token = Materials.readToken(fID, matFileName);
          if isempty(strtrim(token))
            token = Materials.readToken(fID, matFileName);
            assert(~isempty(token),'Multiple newlines after property %s\nIn material file %s', propName, matFileName);
            propName = token;
          elseif strcmp(token,'End')
            break
          else
            error('Wrong number of rows in property %s\nIn material file %s', propName, matFileName)
          end
%           if strcmp(token,'End')
%             break
%           end
%           assert(isempty(str),'Wrong number of rows in material %s\nIn file %s', token, matFileName);
      end
      obj.db(matID) = matProp;
      fclose(fID);
    end
      
      
      
%       % flBlock: flag reporting the reading status
%       %          0 -> the material block has been read correctly
%       %          1 -> error while reading the material block
%       %               (end-of-file before the 'End' statement
%       %          2 -> end of file reached
%       flBlock = 0;
%       % Number of blocks counter
%       nBlock = 0;
%       while (flBlock == 0)
%         % Update the counter
%         nBlock = nBlock + 1;
%         % Reading the material block
%         [flBlock,block] = Materials.readBlock(fid);
%         if flBlock == 0
%           if ~isnan(str2double(block(1)))
%               error('The first entry in block %d of Materials file %s must be strings',nBlock,matFileName);
%           end
%           % Calling the specific material class based on the
%           % material name
%           switch lower(block(1))
%             case 'elastic'
%               obj.db(nBlock) = Elastic(block);
%             case 'hypoelastic'
%               obj.db(nBlock) = HypoElastic(block);
%             case 'transvelastic'
%               obj.db(nBlock) = TransvElastic(block);
%             case 'porousrock'
%               obj.db(nBlock) = PorousRock(block);
%             case 'fluid'
%               obj.db(nBlock) = Fluid(block);
%             otherwise
%               error('Material %s not available (block %d of Materials)', block(1),nBlock);
%           end
%         elseif flBlock == 1
%           error('Error encountered while reading block %d of Material file',nBlock);
%         end
%       end
%       fclose(fid);
%     end
    
    % Reading boundary input file
    function readInputFiles(obj, fListName)
      fID = fopen(fListName,'r');
      Materials.checkFileOpening(fID, fListName);
      matFileName = Materials.readToken(fID, fListName);
      matID = 0;
      while ~strcmp(matFileName,'End')
        matID = matID + 1;
        readInputFile(obj, matFileName, matID);
        matFileName = Materials.readToken(fID, fListName);
      end
      fclose(fID);
    end
  end

  methods (Static = true)
    function checkFileOpening(fID, fName)
      if fID == -1
        error('File %s does not exist in the directory',fName);
      end
    end
    
    % Read the next token and check for eof
    function token = readToken(fID, fName)
      if feof(fID)
        error('Unexpected end of file %s before End statement',fName);
      end
      token = sscanf(fgetl(fID), '%s', 1);
    end
    
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