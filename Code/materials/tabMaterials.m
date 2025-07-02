classdef tabMaterials < handle
  % General class for materials assigned with tabular input format
  % Only Elastic class now support a tabular input format

  properties (Access = private)
    solid
    fluid
  end

  methods (Access = public)
    % Class constructor method   
    function obj = tabMaterials(model,mesh,fName)
      obj.readInputFile(model,mesh,fName)
    end

    % Get the material defined by matIdentifier and check if it is a
    % key of the Map db
    function mat = getMaterial(obj,~)
       mat = obj.solid;
      % 
    end

    function mat = getFluid(obj)
        mat = obj.fluid;
    end
    
    function varargout = computeSwAnddSw(obj,mesh,pkpt)
      % varargout{1} -> Sw
      % varargout{2} -> dSw
      % varargout{2} -> d2Sw
      % if 5<1
      varargout{1} = zeros(mesh.nCells,1);
      if nargout == 2
        varargout{2} = zeros(mesh.nCells,1);
      end
      if nargout == 3
          varargout{2} = zeros(mesh.nCells,1);
          varargout{3} = zeros(mesh.nCells,1);
      end
      for m = 1:mesh.nCellTag
          isElMat = mesh.cellTag == m;
          p = pkpt(isElMat);
          if nargout == 1
              varargout{1}(isElMat) = obj.getMaterial(m).CapillaryCurve.interpTable(p);
          elseif nargout == 2
              [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
          elseif nargout == 3
              [varargout{1}(isElMat), varargout{2}(isElMat), varargout{3}(isElMat)] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
          end
        Swr = obj.getMaterial(m).PorousRock.getWaterResSat();
        varargout{1}(isElMat) = Swr + (1-Swr)*varargout{1}(isElMat);
        if nargout > 1
            varargout{2}(isElMat) = (1-Swr)*varargout{2}(isElMat);
        end
        if nargout > 2
            varargout{3}(isElMat) = (1-Swr)*varargout{3}(isElMat);
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
      mu = obj.getFluid().getDynViscosity();
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

    function [D, sigma, status] = updateMaterial(obj, cTag, sigma, epsilon, dt, status, el, t)
       mat = getMaterial(obj,cTag).ConstLaw;
       [D, sigma, status] = mat.getStiffnessMatrix(sigma, epsilon, dt, status, el);
    end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, model, mesh, matFileName)
      fID = Materials.openReadOnlyFile(matFileName);
      %
      assert(~feof(fID),'No material properties have been assigned in %s',matFileName);
      propName = readToken(fID, matFileName);
      while true
        % Calling the specific material class based on the
        % material name
        switch propName
          case 'Elastic'
            matProp.ConstLaw = Elastic(fID, matFileName, mesh);
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
            obj.fluid = Fluid(fID, matFileName);
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
      obj.solid = matProp;
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