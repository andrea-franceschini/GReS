classdef Materials < handle
  % MATERIAL - General material class
  % 
  % Details - Store capabilities:
  % - for solid materials, you can have only one type material for each
  % subdomain, but you can have the same material defined in more than one
  % subdomain.
  % - for fluid materials (phase), you can have one material define in
  % diferent subdomain, and the subdomain can have more than one material
  % (limit by then physics defined for that subdomain).
  % 
  % This approx still leave the possibility of two subdomain adjacent to
  % each other have different fluid material defined, which can represents
  % a different physics for this subdomains (One that combine all the
  % different phases present).

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
        % varargout{3} -> d2Sw
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
            Sws = obj.getMaterial(m).PorousRock.getMaxSaturation();
            Swr = obj.getMaterial(m).PorousRock.getResidualSaturation();
            if nargout == 1
                [varargout{1}(isElMat), ~, ~] = obj.getMaterial(m).Curves.computeSwAnddSw(p,Swr,Sws);
                % varargout{1}(isElMat) = obj.getMaterial(m).CapillaryCurve.interpTable(p);
                % varargout{1}(isElMat) = Swr + (Sws-Swr)*varargout{1}(isElMat);
            elseif nargout == 2
                [varargout{1}(isElMat), varargout{2}(isElMat), ~] = obj.getMaterial(m).Curves.computeSwAnddSw(p,Swr,Sws);
                % [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).CapillaryCurve.interpTable(p);
                % varargout{1}(isElMat) = Swr + (Sws-Swr)*varargout{1}(isElMat);
                % varargout{2}(isElMat) = (Sws-Swr)*varargout{2}(isElMat);
            elseif nargout == 3
                [varargout{1}(isElMat), varargout{2}(isElMat), varargout{3}(isElMat)] = obj.getMaterial(m).Curves.computeSwAnddSw(p,Swr,Sws);

                % unitSys = 9.8066e2; % SI-9.8066e3 cm-9.8066e2
                % mx = -1e3; % SI-1e1 cm-1e4
                % mn = -7.5e+1; % SI-7.5e-1 cm-7.5e1
                % xx = linspace(mx*unitSys,mn*unitSys,400);
                % [Sa, dSa, ddSa] = obj.getMaterial(m).Curves.computeSwAnddSw(xx,Swr,Sws);
                % figure();
                % hold on;
                % plot(xx/unitSys,Sa,'b', 'LineWidth', 2, 'MarkerSize', 10);
                % % plot(xx,Sa,'b', 'LineWidth', 2, 'MarkerSize', 10);
                % xlabel('Pressure');
                % ylabel('Saturation');
                % xlim([mx,mn]);
                % ylim([Swr,Sws]);
                % ylim([0.102,0.2]);
                % grid on
            end
        end
    end
    
    function varargout = computeLwAnddLw(obj,mesh,upElem,pkpt)
      % Consider moving this method to a more appropriate location
      % varargout{1} -> lw
      % varargout{2} -> dlw
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
          [varargout{1}(isElMat), ~] = obj.getMaterial(m).Curves.computeRelativePermeability(p);
          % varargout{1}(isElMat) = obj.getMaterial(m).RelativePermCurve.interpTable(p);
        elseif nargout == 2
          [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).Curves.computeRelativePermeability(p);
          % [varargout{1}(isElMat), varargout{2}(isElMat)] = obj.getMaterial(m).RelativePermCurve.interpTable(p);          
          % 
          % unitSys = 9.8066e2; % SI-9.8066e3 cm-9.8066e2
          % mx = -1e3; % SI-1e1 cm-1e4
          % mn = -7.5e+1; % SI-7.5e-1 cm-7.5e1
          % xx = linspace(mx*unitSys,mn*unitSys,400);
          % [ka, dka] = obj.getMaterial(m).Curves.computeRelativePermeability(xx);
          % figure();
          % hold on;
          % plot(xx/unitSys,ka,'b', 'LineWidth', 2, 'MarkerSize', 10);
          % % plot(-xx,dkt,'r', 'LineWidth', 2, 'MarkerSize', 10);
          % xlabel('Pressure');
          % ylabel('Relative Permeability');
          % legend({'Analytical'},'Location','northwest');
          % xlim([mx,mn]);
          % grid on;
        
        end
      end
      % if (max(p)>-75)
      %     % pause;
      % end

      mu = obj.getFluid().getDynViscosity();
      varargout{1} = varargout{1}/mu;
      varargout{2} = varargout{2}/mu;
    end

    function [krel, dkrel] = computeRelativePermeability(obj,pres,mat)
        %COMPUTERELATIVEPERMEABILITY compute the relative permeability of
        % a material by a given pressure.
        [krel, dkrel] = obj.getMaterial(mat).Curves.computeRelativePermeability(pres);
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
          case 'VanGenuchten'
            matProp.Curves = VanGenuchten(fID, matFileName);
          case 'Fluid'
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

        % SOLID - Creating the datbase.
        % obj.solidDB = dictionary(uint64.empty,struct.empty);
        obj.solidDB = containers.Map('KeyType','uint64','ValueType','any');
        % Store the information for the Solid class.
        nsolmat = length(fdata.Solid);
        obj.solidMat2Cell = repmat(struct('matID', [], 'cellTag', []), nsolmat, 1);
        for mat=1:nsolmat
            if isnumeric(fdata.Solid(mat).cellTag)
                cellTags = fdata.Solid(mat).cellTag;
            else
                cellTags = sscanf(fdata.Solid(mat).cellTag, '%i,');
            end
            obj.solidMat2Cell(mat).matID = mat;
            obj.solidMat2Cell(mat).cellTag = cellTags;
            [~, ~, extension] = fileparts(fdata.Solid(mat).path);
            switch extension
                case '.xml'
                    material = Materials.readXMLSolid(fdata.Solid(mat).path);
                otherwise
                    % obj.readInputFiles(model,fListName);
            end
            obj.solidDB(mat) = material;
        end
        % Make the map with with subdomain(cellTag) and it's material.
        % obj.solidCell2Mat = dictionary(uint64.empty,uint64.empty);
        obj.solidCell2Mat = containers.Map('KeyType','uint64','ValueType','uint64');
        for mat=1:nsolmat
            temp=obj.solidMat2Cell(mat).cellTag;
            for cells=1:length(temp)
                if ~isKey(obj.solidCell2Mat,temp(cells))
                    obj.solidCell2Mat(temp(cells))=mat;
                else
                    if obj.solidCell2Mat(temp(cells)) ~= mat
                        error('More than one type of solid for the same subdomain!');
                    end
                end
            end
        end

        % FLUID - Creating the datbase.
        % obj.fluidDB = dictionary(uint64.empty,struct.empty);
        obj.fluidDB = containers.Map('KeyType','uint64','ValueType','any');
        % Store the information for the fluid class.
        nflumat = length(fdata.Fluid);
        obj.fluidMat2Cell = repmat(struct('matID', 1, 'cellTag', 1), nflumat, 1);
        for mat=1:nflumat
            fFlowNeed = model.isSinglePhaseFlow() || model.isVariabSatFlow();
            if fFlowNeed && isfield(fdata,'Fluid')
                if isnumeric(fdata.Fluid(mat).cellTag)
                    cellTags = fdata.Fluid(mat).cellTag;
                else
                    cellTags = sscanf(fdata.Fluid(mat).cellTag, '%i,');
                end
                obj.fluidMat2Cell(mat).matID = mat;
                obj.fluidMat2Cell(mat).cellTag = cellTags;
                switch extension
                    case '.xml'
                        material = Materials.readXMLFluid(fdata.Fluid(mat).path);
                    otherwise
                        % obj.readInputFiles(model,fListName);
                end
                obj.fluidDB(mat) = material;
            end
        end
        % Make the map with with subdomain(cellTag) and it's material.
        obj.fluidCell2Mat = containers.Map('KeyType','uint64','ValueType','any');
        scells = unique(vertcat(obj.fluidMat2Cell.cellTag));
        for cel=1:length(scells)
            mats=[];
            for mat=1:nflumat
                if ismember(scells(cel),obj.fluidMat2Cell(mat).cellTag)
                    mats=[mats mat];
                end
            end
            obj.fluidCell2Mat(scells(cel)) = mats;
        end
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