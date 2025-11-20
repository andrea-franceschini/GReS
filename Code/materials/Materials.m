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
    function obj = Materials(input)
      obj.db = containers.Map('KeyType','double','ValueType','any');
      % Calling the function to read input data from file
      obj.readInputFile(input);
    end

    % Get the material defined by matIdentifier and check if it is a
    % key of the Map db
    function mat = getMaterial(obj,cellTag)
      %
      % The preliminary check whether matID is key of db has been commented
      % since it is highly expensive
      %       if (obj.db.isKey(matID))
      mat = obj.db(obj.matMap(cellTag));
      %       else
      %       Displaying error message if matIdentifier is not a key
      %       of the map
      %         error('Material %s not present', matID);
      %       end
    end

    function fluidMat = getFluid(obj)
      % fluid material is always stored as the last one in the database
      fluidMat = obj.db(max(cell2mat(keys(obj.db))));
    end

    function [status] = initializeStatus(obj,cTag,sigma)
      mat = obj.getMaterial(cTag).ConstLaw;
      [status] = mat.initializeStatus(sigma);
    end

    function [D, sigma, status] = updateMaterial(obj, cTag, sigma, epsilon, dt, status, el, t)
      % constitutive update
      mat = obj.getMaterial(cTag).ConstLaw;
      [D, sigma, status] = mat.getStiffnessMatrix(sigma, epsilon, dt, status, el);
    end

    % Destructor
    % function delete(obj)
    %   remove(obj.db,keys(obj.db));
    % end
  end

  methods (Access = private)
    % Reading the input file by material blocks
    function readInputFile(obj, input)

      if ~isstruct(input)
        input = readstruct(input,AttributeSuffix="");
      end

      if isfield(input,"Materials")
        input = input.Materials;
      end

      if isfield(input,"fileName")
        assert(isscalar(fieldnames(inputstruct)),"FileName, " + ...
          " must be a unique parameter.");
        input = readstruct(input.fileName,AttributeSuffix="");
      end

      if isfield(input,"Fluid")
         fluid = Fluid(input.Fluid);
         fluid.name = getXMLData(input.Fluid,[],"name");
      end

      nSolid = 0;
      if isfield(input,"Solid")
        nSolid = numel(input.Solid);
        cTags = [input.Solid.cellTags];
        if isstring(cTags)
          cTags = str2num(strjoin(cTags));
        end
        maxCellTag = max(cTags);
        obj.matMap = zeros(maxCellTag,1);
        for i = 1:nSolid
          mat.name = getXMLData(input.Solid(i),[],"name");
          cellTags = getXMLData(input.Solid(i),[],"cellTags");
          if any(obj.matMap(cellTags))
            existingCellTags = cellTags(obj.matMap(cellTags)~=0);
            gres_log().error("Multiple materials assigned to cellTags %s",...
              sprintf("%i ", existingCellTags));
          end
          obj.matMap(cellTags) = i;
          if isfield(input.Solid(i),"Constitutive")
            if ~ismissing(input.Solid(i).Constitutive)
              constLaws = fieldnames(input.Solid(i).Constitutive);
              % assumes that the XML field has the same name as the
              % constitutive law class
              mat.ConstLaw = feval(constLaws{1},input.Solid(i).Constitutive);
            end
          end
          if isfield(input.Solid(i),"PorousRock")
            if ~ismissing(input.Solid(i).PorousRock)
              mat.PorousRock = PorousRock(input.Solid(i).PorousRock);
            end
          end
          if isfield(input.Solid(i),"Curves")
            if ~ismissing(input.Solid(i).Curves)
              mat.Curves = VanGenuchten(input.Solid(i).Curves);
              mat.Curves.betaCorrection(getFluidSpecWeight(fluid));
            end
          end
          obj.db(i) = mat;
        end
      end

      if ~all(obj.matMap)
        % check if some cell tags do not have any material assigned
        missingCellTags = find(obj.matMap==0);
        error("Missing material for cellTags %s",...
          sprintf("%i ", missingCellTags));
      end

      if isfield(input,"Fluid")
        obj.db(nSolid+1) = fluid;
      end

    end


  end

end