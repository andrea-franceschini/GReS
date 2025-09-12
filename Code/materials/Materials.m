classdef Materials < handle
  % MATERIAL - General material class

  properties (Access = public)
    % Creation of a Dictionary object (faster than map)
    %db = configureDictionary("double","struct");
    % configureDictionary is not supported before 2023b
    db
    matMap
  end

  methods (Access = public)
    % Class constructor method
    function obj = Materials(model,fListName)
      obj.db = containers.Map('KeyType','double','ValueType','any');
      % Calling the function to read input data from file
      obj.matMap = zeros(100,100);
      [~, ~, extension] = fileparts(fListName);
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

    function [status] = initializeStatus(obj,cTag,sigma)
      mat = obj.getMaterial(cTag).ConstLaw;
      [status] = mat.initializeStatus(sigma);
    end

    function [D, sigma, status] = updateMaterial(obj, cTag, sigma, epsilon, dt, status, el, t)
      mat = obj.getMaterial(cTag).ConstLaw;
      [D, sigma, status] = mat.getStiffnessMatrix(sigma, epsilon, dt, status, el);
    end

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
          case 'VanGenuchten'
            matProp.Curves = VanGenuchten(fID, matFileName);
          case 'Fluid'
            assert(isempty(cellTags),'Fluid material do not require cellTag')
            matProp = Fluid(fID, matFileName);
          otherwise
            error('Error in file %s\nMaterial property %s not available', matFileName, propName);
        end
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

      % Test necessity of beta correction in Van Genuchten object.
      if model.isVariabSatFlow()
        if obj.getFluid().isvalid || obj.getMaterial(1).Curves.isvalid
          if obj.getMaterial(1).Curves.needBetaCor()
            specWeight = obj.getFluid().getFluidSpecWeight();
            obj.db(1).Curves.betaCorrection(specWeight);
          end
        end
      end
    end

    function readXMLList(obj, model, fListName)
      %READXMLFILE - function to read the material information file in
      %xml and construct the class object.

      fdata = readstruct(fListName,AttributeSuffix="");
      if ~isfield(fdata,'Solid')
        error('Necessary pass material proprieties to at least one type of Solid Material.')
      end

      % Store the information for the Solid class.
      nsolmat = length(fdata.Solid);
      obj.solidMat = repmat(struct('matID', [], 'cellTag', []), nsolmat, 1);
      for mat=1:nsolmat
        if isnumeric(fdata.Solid(mat).cellTag)
          cellTags = fdata.Solid(mat).cellTag;
        else
          cellTags = sscanf(fdata.Solid(mat).cellTag, '%i,');
        end
        obj.solidMat(mat).matID = mat;
        obj.solidMat(mat).cellTag = cellTags;
        readXMLSolid(obj,fdata.Solid(mat).path,mat,cellTags);
      end

      % Make the map with with subdomain(cellTag) and it's material.
      obj.cellsolid = zeros(max(vertcat(obj.solidMat.cellTag)),1);
      for mat=1:nsolmat
        temp=obj.solidMat(mat).cellTag;
        for cells=1:length(temp)
          if obj.cellsolid(mat)==0
            obj.cellsolid(mat)=temp(cells);
          else
            if obj.cellsolid(mat)~=temp(cells)
              error('More than one type of element at the same subdomain!');
            end
          end
        end
      end

      % Check if fluid information is given and the necessity in the physics.
      obj.fluidMat = repmat(struct('matID', 1, 'cellTag', 1), 1, 1);
      fFlowNeed = model.isSinglePhaseFlow() || model.isVariabSatFlow();
      if fFlowNeed && isfield(fdata,'Fluid')
        if isnumeric(fdata.Solid(mat).cellTag)
          cellTags = fdata.Fluid(1).cellTag;
        else
          cellTags = sscanf(fdata.Fluid(1).cellTag, '%i,');
        end
        readXMLFluid(obj,fdata.Fluid(1).path,mat,cellTags)
      end

      % Make the map with with subdomain(cellTag) and it's material.
      obj.cellfluid = zeros(max(vertcat(obj.fluidMat.cellTag)),1);
      for mat=1:nsolmat
        temp=obj.fluidMat(mat).cellTag;
        for cells=1:length(temp)
          if obj.cellfluid(mat)==0
            obj.cellfluid(mat)=temp(cells);
          else
            if obj.cellfluid(mat)~=temp(cells)
              error('More than one type of element at the same subdomain!');
            end
          end
        end
      end
    end

    function readXMLSolid(obj,fileName,material,celltags)
      %readXMLSolid - function to read the material information relate to
      % the solid part of the model.
      fdata = readstruct(fileName,AttributeSuffix="");

      assert(isscalar(matID),['Zero or Multiple materials assigned' ...
        ' to elements with cellTags %i'], cellID)
      isfield(fdata,'Elastic')
    end

    function readXMLFluid(obj,fileName,material,celltags)
      %readXMLFluid - function to read the material information relate to
      % the fluid part of the model.
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

  end
end