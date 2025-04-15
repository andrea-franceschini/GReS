classdef Materials < handle
  % MATERIAL - General material class

  properties (Access = public)
    % Creation of a Dictionary object (faster than map)
    %db = configureDictionary("double","struct"); 
    % configureDictionary is not supported before 2023b
    db 
    matMap
    % *DB - the database for store the proprieties of the material.
    % *Mat2Cell - map the associate the material and which subdomain is defined.
    % *Cell2Mat - map the associate the subdomain and which material is defined.
    solidDB
    solidMat2Cell
    solidCell2Mat
    fluidDB
    fluidMat2Cell
    fluidCell2Mat
  end

  methods (Access = public)
    % Class constructor method   
    function obj = Materials(model,fListName)
      obj.db = containers.Map('KeyType','double','ValueType','any');
      % Calling the function to read input data from file
      obj.matMap = zeros(100,100);
      [~, ~, extension] = fileparts(fListName);
      assert(isempty(extension),['the %s need to have an extension to', ...
                ' be read the material class'], fListName);
      switch extension
          case '.xml'
              obj.readXMLList(model,fListName);
          otherwise
              obj.readInputFiles(model,fListName);
      end
    end

    % Get the material defined by matIdentifier and check if it is a
    % key of the Map db
    function mat = getMaterial(obj,cellID)
      %
      % The preliminary check whether matID is key of db has been commented
      % since it is highly expensive
%       if (obj.db.isKey(matID))
        [matID,~] = find(obj.matMap == cellID);
        assert(isscalar(matID),['Zero or Multiple materials assigned to elements',...
            ' with cellTags %i'], cellID)
        mat = obj.db(matID);
%       else
%       Displaying error message if matIdentifier is not a key 
%       of the map
%         error('Material %s not present', matID);
%       end
    end

    function fluidMat = getFluid(obj)
        %fluid materials corresponds to null rows in materials map. 
        %Fluid are not assigned to any cellTag
        f = find(sum(obj.matMap,2)==0);
        fluidMat = obj.db(f);
    end
    
    function varargout = computeSwAnddSw(obj,mesh,pkpt)
      % Consider moving this method to a more appropriate location
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
      % Consider moving this method to a more appropriate location
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
       mat = obj.getMaterial(cTag).ConstLaw;
       [D, sigma, status] = mat.getStiffnessMatrix(sigma, epsilon, dt, status, el);
    end

    %
    % Destructor
    function delete(obj)
      remove(obj.db,keys(obj.db));
    end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, model, matFileName, matID, cellTags)
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
      if isfield(matProp,'ConstLaw')
          %check if non-fluid material has cellTags assigned
          assert(~isempty(cellTags),"Undefined cellTags for material #%i", matID)
      else
          assert(isempty(cellTags),"Fluid material do not require any cellTag", matID)
      end
      obj.db(matID) = matProp;
      fclose(fID);
    end
         
    % Reading boundary input file
    function readInputFiles(obj, model, fListName)
      fID = Materials.openReadOnlyFile(fListName);
      [cellTags,matFileName] = readTokenList(fID, fListName);
      matID = 0;
      while ~strcmp(matFileName,'End')
        matID = matID + 1;
        readInputFile(obj, model, matFileName, matID, cellTags);
        obj.matMap(matID,1:length(cellTags)) = cellTags;
        [cellTags,matFileName] = readTokenList(fID, fListName);
      end
      obj.matMap = obj.matMap(1:matID,sum(obj.matMap,1)~=0);
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

    function material = readXMLSolid(fileName)
        %readXMLSolid - function to read the material information relate to
        % the solid part of the model.

        material = struct();
        fdata = readstruct(fileName,AttributeSuffix="");
        fnames = fieldnames(fdata);
        for mat=1:length(fnames)
            switch string(fnames(mat))
                case 'Elastic'
                    material.ConstLaw = 0;
                    youngMod = checkField(fdata.Elastic.young,1);
                    nu = fdata.Elastic.poisson;
                    % material.ConstLaw = Elastic(fID, matFileName);
                case 'HypoElastic'
                    % material.ConstLaw = HypoElastic(fID, matFileName);
                case 'TransvElastic'
                    % material.ConstLaw = TransvElastic(fID, matFileName);
                case 'SSCM'
                    % material.ConstLaw = SSCM(fID, matFileName);
                case 'PorousRock'
                    % material.PorousRock = PorousRock(fID, model, matFileName);
                case 'CapillaryCurve'
                    % material.CapillaryCurve = TabularCurve(fID, matFileName);
                case 'RelativePermCurve'
                    % material.RelativePermCurve = TabularCurve(fID, matFileName);
                case 'VanGenuchten'
                    % material.Curves = VanGenuchten(fID, matFileName);
            end
        end
    end



    function material = readXMLFluid(fileName)
        %readXMLFluid - function to read the material information relate to
        % the fluid part of the model.
        material = struct();
    end

    function vec = checkField(stData,npos)
        %CHECKFIELD - function to check if the field inside a struct is a
        % number or a path and return the information.
        if isnumeric(stData)
            vec = stData;
        else
            [~, ~, extension] = fileparts(stData);
            assert(isempty(extension),['the %s need to have an extension to', ...
                ' be read as a list of values. Error in the definition of ', ...
                'XML file for the material class'], stData);
            % read the file.
            vec = 0.;
        end
    end


  end
end