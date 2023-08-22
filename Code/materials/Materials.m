classdef Materials < handle
  % MATERIAL - General material class

  properties (Access = private)
    % Creation of a Map object 
    db = containers.Map('KeyType','double','ValueType','any');
  end

  methods (Access = public)
    % Class constructor method   
    function obj = Materials(model,fListName)
      % Calling the function to read input data from file
      obj.readInputFiles(model,fListName)
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
    
    function varargout = computeSwAnddSw(obj,mesh,pkpt)
      % varargout{1} -> Sw
      % varargout{2} -> dSw
      % if 5<1
      varargout{1} = zeros(mesh.nCells,1);
      if nargout == 2
        varargout{2} = zeros(mesh.nCells,1);
      end
      for m = 1:mesh.nCellTag
        isElMat = mesh.cellTag == m;
        p = pkpt(isElMat);
        if nargout == 1
          varargout{1}(isElMat) = obj.getMaterial(m).CapillaryCurve.interpTable(p);
        elseif nargout == 2
          [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
        end
        Swr = obj.getMaterial(m).PorousRock.getWaterResSat();
        varargout{1}(isElMat) = Swr + (1-Swr)*varargout{1}(isElMat);
        if nargout == 2
          varargout{2}(isElMat) = (1-Swr)*varargout{2}(isElMat);
        end
      end
      % end
      %
      %
      %
      if 5<1
      varargout{1} = ones(mesh.nCells,1);
      if nargout == 2
        varargout{2} = zeros(mesh.nCells,1);
      end
      n = 3.1769;
      pEntry = 2.7840;
      m = 1 - 1/n;
      pkpt = -pkpt;
      isPos = pkpt >= 0;
      SeFun = @(p) (1 + (p./pEntry).^n).^(-m);
      dSeFun = @(p) -m.*(1 + (p./pEntry).^n).^(-m-1).*n./pEntry.*(p./pEntry).^(n-1);
      varargout{1}(isPos) = SeFun(pkpt(isPos));
      Swr = obj.getMaterial(1).PorousRock.getWaterResSat();
      varargout{1} = Swr + (1-Swr)*varargout{1};
      if nargout == 2
        varargout{2}(isPos) = dSeFun(pkpt(isPos));
        varargout{2} = (1-Swr)*varargout{2};
      end
      end
    end
    
    function varargout = computeLwAnddLw(obj,mesh,upElem,pkpt)
      % varargout{1} -> lw
      % varargout{2} -> dlw
      % if 5<1
      nIntFaces = length(upElem);
      varargout{1} = zeros(nIntFaces,1);
      if nargout == 2
        varargout{2} = zeros(nIntFaces,1);
      end
      matUpElem = mesh.cellTag(upElem);
      for m = 1:mesh.nCellTag
        isElMat = matUpElem == m;
        p = pkpt(upElem(isElMat));
        if nargout == 1
          varargout{1}(isElMat) = obj.getMaterial(m).RelativePermCurve.interpTable(p);
        elseif nargout == 2
          [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).RelativePermCurve.interpTable(p);
        end
      end
      mu = obj.getMaterial(mesh.nCellTag+1).getDynViscosity();
      varargout{1} = varargout{1}/mu;
      varargout{2} = varargout{2}/mu;
      % end
      %
      %
      %
      if 5<1
      nIntFaces = length(upElem);
      varargout{1} = ones(nIntFaces,1);
      if nargout == 2
        varargout{2} = zeros(nIntFaces,1);
      end
      n = 3.1769;
      pEntry = 2.7840;
      m = 1 - 1/n;
      pkpt = -pkpt;
      pres = pkpt(upElem);
      isPos = pres >= 0;
      krFun = @(p) (1 + (p/pEntry).^n).^(-2.5.*m) .* ((1 + (p/pEntry).^n).^m - ...
        ((p/pEntry).^n).^m).^2;
      dkrFun = @(p) (-2.5.*m.*((1+(p/pEntry).^n).^(-2.5*m-1)).* ...
          ((1+(p/pEntry).^n).^m - ((p/pEntry).^n).^m).^2 + 2.* ...
          ((1+(p/pEntry).^n).^m - ((p/pEntry).^n).^m).* ...
          (m.*((1+(p/pEntry).^n).^(m-1))-m.*(((p/pEntry).^n).^(m-1))).* ...
          ((1+(p/pEntry).^n).^(-2.5.*m))).*n./pEntry.*(p./pEntry).^(n-1);
      varargout{1}(isPos) = krFun(pres(isPos));
      if nargout == 2
        varargout{2}(isPos) = dkrFun(pres(isPos));
      end
      %
      mu = obj.getMaterial(mesh.nCellTag+1).getDynViscosity();
      varargout{1} = varargout{1}/mu;
      varargout{2} = varargout{2}/mu;
      end
    end
    
    function [status] = initializeStatus(obj,cTag,sigma)
      mat = obj.getMaterial(cTag).ConstLaw;
      [status] = mat.initializeStatus(sigma);
    end
    
    function [D, sigma, status] = updateMaterial(obj, cTag, sigma, epsilon, dt, status, t)
      mat = obj.getMaterial(cTag).ConstLaw;
      [D, sigma, status] = mat.getStiffnessMatrix(sigma, epsilon, dt, status, t);
    end
    
%     function [Sw,dSw,lw,dlw] = computeSwAndLambda(obj,mesh,upElem,pkpt)
%       % kr, dkr -> computed for every internal face
%       % Sw, dSw -> computed for every element
%       nIntFaces = length(upElem);
%       Sw = zeros(mesh.nCells,1);
%       dSw = zeros(mesh.nCells,1);
%       lw = zeros(nIntFaces,1);
%       dlw = zeros(nIntFaces,1);
%       matUpElem = mesh.cellTag(upElem);
%       for m = 1:mesh.nCellTag
%         isElMat = matUpElem == m;
%         p = pkpt(upElem(isElMat));
%         [lw(isElMat), dlw(isElMat)] = obj.getMaterial(m).RelativePermCurve.interpTable(p);
%         clear isElMat p
%         isElMat = mesh.cellTag == m;
%         p = pkpt(isElMat);
%         [Sw(isElMat), dSw(isElMat)] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
%         Swr = obj.getMaterial(m).PorousRock.getWaterResSat();
%         Sw(isElMat) = Swr + (1-Swr)*Sw(isElMat);
%         dSw(isElMat) = (1-Swr)*dSw(isElMat);
%       end
%       mu = obj.getMaterial(obj.mesh.nCellTag+1).getDynViscosity();
%       lw = lw/mu;
%       dlw = dlw/mu;
%     end
    %
    % Destructor
    function delete(obj)
      remove(obj.db,keys(obj.db));
    end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, model, matFileName, matID)
      fID = Materials.openReadOnlyFile(matFileName);
      %
      assert(~feof(fID),'No material properties have been assigned in %s',matFileName);
      propName = readToken(fID, matFileName);
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
          case 'SSCM'
            matProp.ConstLaw = SSCM(fID, matFileName);
          case 'PorousRock'
            matProp.PorousRock = PorousRock(fID, model, matFileName);
          case 'CapillaryCurve'
            matProp.CapillaryCurve = TabularCurve(fID, matFileName);
          case 'RelativePermCurve'
            matProp.RelativePermCurve = TabularCurve(fID, matFileName);
          case 'Fluid'
            matProp = Fluid(fID, matFileName);
          otherwise
            error('Error in file %s\nMaterial property %s not available', matFileName, propName);
        end
%           [str, flEoF] = Materials.readToken(fID);
        token = readToken(fID, matFileName);
        if isempty(strtrim(token))
          token = readToken(fID, matFileName);
          assert(~isempty(token),'Multiple newlines after property %s\nIn material file %s', propName, matFileName);
          propName = token;
        elseif strcmp(token,'End')
          break
        end
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
    function readInputFiles(obj, model, fListName)
      fID = Materials.openReadOnlyFile(fListName);
      matFileName = readToken(fID, fListName);
      matID = 0;
      while ~strcmp(matFileName,'End')
        matID = matID + 1;
        readInputFile(obj, model, matFileName, matID);
        matFileName = readToken(fID, fListName);
      end
      fclose(fID);
    end
  end

  methods (Static = true)
    % Open a read-only file and check  
    function fID = openReadOnlyFile(fName)
      fID = fopen(fName, 'r');
      if fID == -1
        error('File %s does not exist in the directory',fName);
      end
    end

    % Read the material block
%     function [flBlock,block] = readBlock(fid)
%       block = [];
%       % flBegBlock: flag marking the beginning of a material block
%       flBegBlock = false;
%       [flEof,line] = Materials.readLine(fid);
%       % Dealing with the initial blank lines, if any
%       % flEof: end-of-file flag
%       %        1 -> end of file reached
%       %        0 -> otherwise
%       while (isempty(strtrim(line)) && flEof == 0)
%         [flEof,line] = Materials.readLine(fid);
%       end
%       %
%       if flEof == 0
%         flBegBlock = true;
%       end
%       % Reading the material block
%       while (~strcmp(line,'End') && flEof == 0)
%         line = string(strtrim(strtok(line,'%')));
%         block = [block; line];
%         [flEof,line] = Materials.readLine(fid);
%       end
%       %
%       if flEof == 0 && flBegBlock
%         flBlock = 0;   % block read correctly
%       elseif flEof == 1 && flBegBlock
%         flBlock = 1;   % end of file while reading
%       elseif flEof == 1 && ~flBegBlock
%         flBlock = 2;   % end of file reached
%       end
%     end

    % Read the next line and check for eof
%     function [flEof,line] = readLine(fid)
%       flEof = feof(fid);   % end-of-file flag
%       if flEof == 1
%         line = '';
%       else
%         line = strtrim(fgetl(fid));
%       end  
%     end
  end
end